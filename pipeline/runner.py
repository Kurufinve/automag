"""
automag.pipeline.runner
=======================

Executes the stage *instances* of a pipeline.

The runner does not import the stage scripts.  It materialises the work
directories, writes a per-instance ``input.py`` and invokes the existing,
well-tested stage scripts as subprocesses with

* ``cwd``                = ``<workdir>[/<stage>[/<instance>]]``
* ``AUTOMAG_PATH``       = ``<workdir>``          (legacy data root)
* ``AUTOMAG_WORKDIR``    = ``<workdir>``          (new data root)
* ``AUTOMAG_CODE_PATH``  = repository root        (where the code lives)
* ``PYTHONPATH``         = repository root + inherited value

An **instance** is one run of one stage with one calculator
(:class:`~automag.pipeline.config.StageInstance`).  A stage may appear several
times, e.g. ``2_coll`` once with a fast model and once with DFT.  Instances are
scheduled in declaration order.

In the default mode the requested phases are executed in order.  With
``wait=True`` the runner builds a dependency graph from
:mod:`automag.pipeline.stages`, submits a batch of calculations, polls until
they finish (:mod:`automag.pipeline.readiness`) and only then runs the
dependent phases.  This makes ``run_automag.py --wait`` a single unattended
command for the whole pipeline.

Three policies make a long mixed pipeline practical:

``reuse``
    skip an instance whose output already exists (warn instead of recompute);
``on_missing``
    ``"skip"`` skips an instance whose prerequisite results are absent,
    ``"error"`` (default) aborts;
``on_unsupported``
    ``"skip"`` (default for dict instances) skips an instance whose calculator
    cannot perform the stage, e.g. MAE with a spin-lattice model.

.. codeauthor:: Automag contributors
"""
from __future__ import annotations


import os
import pprint
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pipeline.readiness as readiness
from calculators import get_calculator
from common.paths import code_root
from pipeline import execution
from pipeline.config import ConfigError
from pipeline.stages import STAGES, STAGE_ORDER

#: Blocks appended to input.py by scripts (e.g. coupling constants) must
#: survive regeneration of the file.
_APPENDED_BLOCK = re.compile(
    r"(# LINE ADDED BY THE SCRIPT[^\n]*\n[^\n]*)", re.MULTILINE
)

#: Capability flag a calculator interface must set to run a stage.
CAPABILITIES = {
    "lin_response": "supports_ldau",
    "relax": "supports_relaxation",
    "mae": "supports_noncollinear",
}

#: Stages that read the results of a ``2_coll`` instance.
COLL_DEPENDENT = ("relax", "monte_carlo", "mae")

#: Settings that control the runner and must not end up in ``input.py``.
IGNORED_INPUT_KEYS = {"steps", "enabled", "description"}


@dataclass
class Task:
    """One phase of one stage instance."""

    instance: object
    phase: str
    scripts: tuple
    depends: tuple = ()
    wait: bool = False
    wait_before: str | None = None

    @property
    def stage(self):
        return self.instance.stage

    @property
    def key(self):
        return (self.instance.id, self.phase)

    @property
    def label(self):
        return f"{self.instance.id}.{self.phase}"


class PipelineRunner:
    def __init__(self, config, steps=None, stages=None, instances=None,
                 dry_run=False, wait=False, wait_timeout=None,
                 wait_interval=None, wait_stable_seconds=None,
                 wait_for_mc_output=None, mode=None, execution=None,
                 fresh=False, state_file=None, poll_interval=None):
        self.config = config
        self.steps = steps
        self.stages = stages
        self.instances = instances
        self.dry_run = dry_run
        #: ``steps`` (default), ``wait``, ``once`` or ``daemon``.
        self.mode = mode or ("wait" if wait else "steps")
        #: Force one execution mode for every step (``--execution``).
        self.execution_override = execution
        self.fresh = fresh
        self.state_file = state_file
        self.poll_interval = (poll_interval if poll_interval is not None
                              else config.poll_interval)
        self.wait = self.mode in ("wait", "once", "daemon")
        self.wait_timeout = (config.wait_timeout if wait_timeout is None
                             else wait_timeout)
        self.wait_interval = (config.wait_interval if wait_interval is None
                              else wait_interval)
        self.wait_stable_seconds = (config.wait_stable_seconds
                                    if wait_stable_seconds is None
                                    else wait_stable_seconds)
        self.wait_for_mc_output = (config.wait_for_monte_carlo_output
                                   if wait_for_mc_output is None
                                   else wait_for_mc_output)
        self._formula = None
        self._input_written = set()
        self._skips = {}
        self._blocking = set()
        self._hosts = None
        self._validate()

    # ------------------------------------------------------------------
    # Instance selection
    # ------------------------------------------------------------------
    def selected_instances(self):
        """Instances to execute in this invocation, in declaration order."""
        result = list(self.config.instances)
        if self.stages:
            unknown = [name for name in self.stages if name not in STAGE_ORDER]
            if unknown:
                raise ConfigError(
                    f"Unknown stage(s) {unknown}. Known stages: "
                    f"{', '.join(STAGE_ORDER)}.")
            result = [item for item in result if item.stage in self.stages]
        if self.instances:
            known = {item.id for item in self.config.instances}
            unknown = [name for name in self.instances if name not in known]
            if unknown:
                raise ConfigError(f"Unknown instance(s) {unknown}.")
            result = [item for item in result if item.id in self.instances]
        return result

    def formula(self):
        if self._formula is None:
            from pymatgen.core.structure import Structure

            structure = Structure.from_file(str(self.config.structure))
            self._formula = structure.formula.replace(" ", "")
        return self._formula

    def _validate(self):
        selected = self.selected_instances()
        if not selected:
            raise ConfigError("No stage instance selected for this run.")
        for inst in selected:
            stage = STAGES[inst.stage]
            missing = [key for key in stage.required
                       if key not in inst.config]
            if missing:
                raise ConfigError(
                    f"Instance '{inst.id}' ({inst.stage}) is missing required "
                    f"setting(s): {', '.join(missing)}.")

    # ------------------------------------------------------------------
    # Preparation
    # ------------------------------------------------------------------
    def prepare(self):
        workdir = self.config.workdir
        (workdir / "geometries").mkdir(parents=True, exist_ok=True)
        (workdir / "CalcFold").mkdir(parents=True, exist_ok=True)

        destination = workdir / "geometries" / self.config.structure.name
        source = self.config.structure
        if source.resolve() != destination.resolve():
            shutil.copy2(source, destination)

        for inst in self.selected_instances():
            inst.work_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Input file generation
    # ------------------------------------------------------------------
    def write_stage_input(self, inst):
        inst.work_dir.mkdir(parents=True, exist_ok=True)
        lines = [
            "# Auto-generated by automag.pipeline from",
            f"# {self.config.path}",
            f"# instance '{inst.id}': stage {inst.stage}, "
            f"calculator {inst.calculator}",
            "# Edit the controlling file, not this one.",
            "",
        ]
        for key, value in inst.config.items():
            if key in IGNORED_INPUT_KEYS or key.startswith("_"):
                continue
            rendered = pprint.pformat(value, width=100, sort_dicts=False)
            lines.append(f"{key} = {rendered}")
            lines.append("")

        input_file = inst.work_dir / "input.py"
        appended = ""
        if input_file.is_file():
            appended = "\n".join(_APPENDED_BLOCK.findall(input_file.read_text()))
        if appended:
            lines.append(appended)
            lines.append("")

        input_file.write_text("\n".join(lines))
        return input_file

    # ------------------------------------------------------------------
    # Environment
    # ------------------------------------------------------------------
    def environment(self):
        env = os.environ.copy()
        root = str(code_root())
        workdir = str(self.config.workdir)
        env["AUTOMAG_CODE_PATH"] = root
        env["AUTOMAG_WORKDIR"] = workdir
        # Keep the historical variable pointing at the data root so that the
        # existing stage scripts continue to work unchanged.
        env["AUTOMAG_PATH"] = workdir
        pythonpath = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = root + (os.pathsep + pythonpath if pythonpath else "")
        if self.config.calculator_plugins:
            env["AUTOMAG_CALCULATORS"] = ",".join(self.config.calculator_plugins)
        # Stage scripts hand this to SubmitManual: 0 = write inputs only.
        env["AUTOMAG_AUTO_SUBMIT"] = (
            "0" if self.config.common.get("auto_submit") is False else "1")
        return env

    # ------------------------------------------------------------------
    # Policies: unsupported calculator, reuse, missing prerequisites
    # ------------------------------------------------------------------
    def _supported(self, inst):
        attribute = CAPABILITIES.get(inst.stage)
        if attribute is None:
            return True
        try:
            interface = get_calculator(
                inst.calculator, plugins=self.config.calculator_plugins)
        except Exception as error:  # noqa: BLE001 - report a clean message
            raise ConfigError(
                f"Cannot load calculator '{inst.calculator}' for instance "
                f"'{inst.id}': {error}") from error
        return bool(getattr(interface, attribute, False))

    def _stage_done(self, inst):
        """Return True when the output of ``inst`` already exists on disk."""
        from common.results import most_stable_state

        workdir = self.config.workdir
        formula = self.formula()
        suffix = self.config.struct_suffix
        calc = inst.calculator

        if inst.stage == "coll":
            try:
                most_stable_state(workdir, formula, suffix, calc)
                return True
            except FileNotFoundError:
                return False
        if inst.stage == "relax":
            return bool(list(inst.work_dir.glob("relaxed_*.json")))
        if inst.stage == "monte_carlo":
            output = inst.work_dir / "output"
            return output.is_file() and output.stat().st_size > 0
        if inst.stage == "mae":
            return bool(list(inst.work_dir.glob("MAE_curve_*.txt")))
        if inst.stage == "conv_tests":
            mode = inst.config.get("mode", "kgrid")
            path = workdir / "CalcFold" / f"{formula}{suffix}_{mode}_{calc}.txt"
            return path.is_file()
        if inst.stage == "lin_response":
            path = workdir / "CalcFold" / f"{formula}{suffix}_charges_{calc}.txt"
            return path.is_file()
        return False

    def _coll_data_available(self, inst):
        """True when 2_coll results for ``inst.coll_calculator`` exist."""
        from common.results import most_stable_state

        try:
            most_stable_state(self.config.workdir, self.formula(),
                              self.config.struct_suffix, inst.coll_calculator)
            return True
        except FileNotFoundError:
            return False

    def plan_skips(self, instances):
        """Decide which instances must be skipped, and why."""
        by_id = {item.id: item for item in self.config.instances}
        selected = {item.id for item in instances}
        skips, blocking = {}, set()

        for inst in instances:
            reason, policy, blocks = None, None, False
            if not self._supported(inst):
                reason = (f"calculator '{inst.calculator}' cannot perform the "
                          f"'{inst.stage}' stage")
                policy = inst.on_unsupported
                blocks = True
            elif inst.reuse and self._stage_done(inst):
                reason = "its output already exists and 'reuse' is enabled"
            else:
                for dep_id in inst.depends_on:
                    if dep_id in blocking:
                        reason = f"prerequisite instance '{dep_id}' was skipped"
                        policy = inst.on_missing
                        blocks = True
                        break
                    dep = by_id[dep_id]
                    if dep_id in selected or self._stage_done(dep):
                        continue
                    reason = (f"the results of prerequisite instance "
                              f"'{dep_id}' ({dep.stage}, {dep.calculator}) "
                              f"are missing")
                    policy = inst.on_missing
                    blocks = True
                    break
                if (reason is None and inst.stage in COLL_DEPENDENT
                        and not inst.depends_on
                        and not self._coll_data_available(inst)):
                    # No coll instance is part of this run and none is
                    # recorded as a dependency: fall back to the data on disk.
                    reason = (f"no 2_coll results for calculator "
                              f"'{inst.coll_calculator}'")
                    policy = inst.on_missing
                    blocks = True

            if reason is None:
                continue
            if policy == "error":
                raise ConfigError(
                    f"Cannot run instance '{inst.id}' ({inst.stage}, "
                    f"{inst.calculator}): {reason}.")
            skips[inst.id] = reason
            if blocks:
                blocking.add(inst.id)
            print(f"WARNING: skipping instance '{inst.id}' "
                  f"({inst.stage}, {inst.calculator}): {reason}.")
            if inst.stage in COLL_DEPENDENT:
                self._warn_alternatives(inst)
        return skips, blocking

    def _warn_alternatives(self, inst):
        """Point the user at collinear results they did not ask to analyse."""
        from common.results import available_coll_calculators

        try:
            found = available_coll_calculators(
                self.config.workdir, self.formula(),
                self.config.struct_suffix)
        except Exception:  # noqa: BLE001 - purely advisory
            return
        if not found:
            return
        print(f"  note: 2_coll results exist for calculator(s): "
              f"{', '.join(sorted(found))}.")
        declared = {item.coll_calculator for item in self.config.instances
                    if item.stage == inst.stage}
        extra = sorted(found - declared)
        if extra:
            print(f"  note: add a '{inst.stage}' instance with "
                  f"'calculator'/'coll_calculator' in {extra} to analyse "
                  f"those results.")

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------
    def run(self):
        self.prepare()
        instances = self.selected_instances()
        self._skips, self._blocking = self.plan_skips(instances)
        runnable = [item for item in instances
                    if item.id not in self._skips]
        if not runnable:
            print("Nothing to do: every instance was skipped.")
            return 0
        if self.mode in ("wait", "once", "daemon"):
            return self.run_scheduled(runnable)
        env = self.environment()
        for inst in runnable:
            self.run_instance(inst, env)
        return 0

    # ------------------------------------------------------------------
    # State driven execution (--wait / --once / --daemon)
    # ------------------------------------------------------------------
    def state_path(self):
        """Location of the SQLite file holding the pipeline state."""
        if self.state_file is not None:
            path = Path(self.state_file).expanduser()
            return path if path.is_absolute() else self.config.workdir / path
        if self.config.state_file is not None:
            return Path(self.config.state_file)
        return self.config.workdir / "automag_state.sqlite"

    def run_scheduled(self, runnable):
        from pipeline.scheduler import Scheduler
        from pipeline.state import StateStore

        path = self.state_path()
        if self.dry_run:
            state = StateStore(":memory:")
        else:
            if self.fresh and Path(path).is_file():
                Path(path).unlink()
            state = StateStore(path)
        self._state = state
        scheduler = Scheduler(
            self, state, mode=self.mode, poll_interval=self.poll_interval,
            fresh=self.fresh,
            execution_override=self.execution_override)
        scheduler.prepare()
        return scheduler.run()

    def run_instance(self, inst, env):
        stage = STAGES[inst.stage]
        selection = inst.config.get("steps", self.steps)
        scripts = stage.resolve_steps(selection)
        inst.work_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n=== Instance '{inst.id}' "
              f"(stage {inst.stage}, calculator {inst.calculator}) ===")
        print(f"  dir: {inst.work_dir}")
        input_file = self.write_stage_input(inst)
        print(f"  input: {input_file}")

        for script in scripts:
            self._run_script(inst, script, env)

    # ------------------------------------------------------------------
    # Dependency-aware (--wait) execution
    # ------------------------------------------------------------------
    def _planned_phases(self, inst):
        stage = STAGES[inst.stage]
        selection = inst.config.get("steps", self.steps)
        full = self.wait and selection is None
        return stage.phases(selection, full=full)

    def build_tasks(self, instances):
        planned = {item.id: self._planned_phases(item) for item in instances}
        tasks = []
        for inst in instances:
            stage = STAGES[inst.stage]
            phases = planned[inst.id]
            for index, phase in enumerate(phases):
                depends = []
                for dep_stage, dep_phase in stage.depends.get(phase, ()):
                    if dep_stage == inst.stage:
                        depends.append((inst.id, dep_phase))
                if index == 0:
                    for dep_id in inst.depends_on:
                        dep_phases = planned.get(dep_id)
                        depends.append((dep_id,
                                        dep_phases[-1] if dep_phases else None))
                wait_before = None
                if inst.stage == "monte_carlo" and phase == "postprocess":
                    wait_before = "mc_output"
                tasks.append(Task(
                    instance=inst,
                    phase=phase,
                    scripts=stage.steps[phase],
                    depends=tuple(depends),
                    wait=phase in stage.wait_after,
                    wait_before=wait_before,
                ))
        return tasks, planned

    def _wait_before(self, task, run_selected=False):
        if task.wait_before != "mc_output":
            return True

        output = task.instance.work_dir / "output"
        if output.is_file() and output.stat().st_size > 0:
            return True
        if self.dry_run:
            print(f"  (dry-run) would check for {output}")
            return True
        if run_selected:
            raise ConfigError(
                f"The monte_carlo run phase of '{task.instance.id}' finished "
                f"but {output} was not produced. Check 'mc_backend', the run "
                "command and the output file name.")
        if not self.wait_for_mc_output:
            print(
                f"WARNING: {output} not found. The Monte Carlo simulation "
                "was not run by the pipeline (the 'run' phase was skipped). "
                "Skipping the monte_carlo postprocess phase; set "
                "'wait_for_monte_carlo_output = True' to wait for an "
                "externally produced output.")
            return False

        readiness.wait_for_file(
            output, timeout=self.wait_timeout, interval=self.wait_interval,
            label=f"Monte Carlo output {output}")
        return True

    def host_for(self, inst):
        """Resolved host of one stage instance."""
        from pipeline import hosts as hosts_module

        if self._hosts is None:
            self._hosts = hosts_module.resolve_hosts(
                self.config, local_workdir=self.config.workdir)
        return hosts_module.host_for(inst, self._hosts)

    def _run_script(self, inst, script, env):
        directory = STAGES[inst.stage].directory
        script_path = code_root() / directory / script
        if not script_path.is_file():
            raise FileNotFoundError(f"Stage script not found: {script_path}")
        host = self.host_for(inst)
        mode = execution.resolve_mode(
            inst, None, host=host, override=self.execution_override)
        print(f"  running: {directory}/{script} ({mode} @{host.name})")
        if self.dry_run:
            return
        if mode == "inprocess":
            code = execution.run_inprocess(script_path, inst.work_dir, env)
            if code != 0:
                raise subprocess.CalledProcessError(code, str(script_path))
            return
        if not host.is_local:
            step_env = {key: value for key, value in env.items()
                        if key.startswith("AUTOMAG_") or key == "PYTHONPATH"}
            code, _, _ = execution.run_sync(
                host, execution.automag_command(host, f"{directory}/{script}"),
                inst.work_dir, env=step_env)
            if code != 0:
                raise subprocess.CalledProcessError(code, script_path.name)
            return
        # The stage scripts do ``from input import *``.  ``python <abs path>``
        # puts the *script's* directory on sys.path, not the cwd, so the
        # instance work directory has to be added explicitly.
        task_env = dict(env)
        parts = [str(inst.work_dir), str(code_root())]
        if env.get("PYTHONPATH"):
            parts.append(env["PYTHONPATH"])
        task_env["PYTHONPATH"] = os.pathsep.join(parts)

        if getattr(sys, "frozen", False):
            # A one-file build unpacks to a fresh temporary directory on every
            # start, so the child cannot be given this process's absolute path.
            command = [sys.executable, "--run-script", f"{directory}/{script}"]
        else:
            command = [sys.executable, str(script_path)]
        subprocess.run(
            command,
            cwd=str(inst.work_dir),
            env=task_env,
            check=True,
        )


def stage_summary():
    """Human-readable description of all stages and their phases."""
    lines = []
    for name in STAGE_ORDER:
        stage = STAGES[name]
        phases = ", ".join(stage.steps)
        lines.append(f"{name:12s} {stage.directory:14s} phases: {phases}")
    return "\n".join(lines)
