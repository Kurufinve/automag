"""
automag.pipeline.scheduler
==========================

State driven execution of the stage instances.

The scheduler turns the pipeline into a resumable state machine stored in
:class:`pipeline.state.StateStore` (SQLite).  Every unit of work is a *step*:
one phase (``submit``, ``postprocess``, ``axes``, ...) of one stage instance.

Two run modes are built on top of it:

``--once`` (one-shot)
    advance the pipeline until a step needs to wait for external programs,
    record the started jobs (host, scheduler id / pid, calculation list) and
    exit.  Running the same command again continues where it stopped.

``--daemon`` / ``--wait`` (continuous)
    keep polling the started jobs and start the next step as soon as the
    previous one is complete.  ``--wait`` exits when the pipeline is done,
    ``--daemon`` can optionally stay alive and watch the input file.

Changed settings are detected per step (``config_hash``): only the steps whose
effective settings changed are reset, every other step - and every individual
calculation that already produced a result - is kept.

.. codeauthor:: Automag contributors
"""
from __future__ import annotations


import hashlib
import json
import shlex
import sys
import time
from pathlib import Path

import pipeline.readiness as readiness
from common.paths import code_root
from pipeline import execution
from pipeline.config import ConfigError
from pipeline.stages import STAGES

#: Settings that do not influence the physics of a step and therefore must not
#: trigger a restart when they change.
HASH_EXCLUDED = {
    "description", "steps", "enabled", "reuse", "on_missing", "on_unsupported",
    "max_parallel", "max_parallel_tasks", "parallel", "max_attempts",
    "execution", "host", "work_dir", "state_file", "poll_interval",
}


def step_key(instance, phase):
    return f"{instance.id}.{phase}"


def config_hash(instance, phase, scripts, host_name=None, structure=None):
    """Stable hash of everything that decides what a step computes."""
    settings = {key: value for key, value in instance.config.items()
                if key not in HASH_EXCLUDED and not key.startswith("_")}
    structure_stamp = None
    if structure is not None:
        try:
            stat = Path(structure).stat()
            structure_stamp = [str(structure), stat.st_size, int(stat.st_mtime)]
        except OSError:
            structure_stamp = str(structure)
    payload = {
        "instance": instance.id,
        "stage": instance.stage,
        "phase": phase,
        "calculator": instance.calculator,
        "coll_calculator": instance.coll_calculator,
        "scripts": list(scripts),
        "host": host_name,
        "settings": settings,
        "structure": structure_stamp,
    }
    text = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha1(text.encode()).hexdigest()[:16]


class Scheduler:
    """Execute the tasks of :meth:`pipeline.runner.PipelineRunner.build_tasks`."""

    def __init__(self, runner, state, mode="once", poll_interval=None,
                 fresh=False, execution_override=None):
        self.runner = runner
        self.state = state
        self.config = runner.config
        self.mode = mode
        self.fresh = fresh
        self.execution_override = execution_override
        self.poll_interval = float(
            poll_interval or self.config.poll_interval or 30)
        self._input_written = set()
        self._skipped = {}
        self._jobs = {}
        self._hashes = {}
        self._tasks = []
        self._planned = {}
        self._structure = None
        self._hosts = None
        self.dry_run = runner.dry_run

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------
    def _host(self, instance):
        from pipeline import hosts as hosts_module

        if self._hosts is None:
            self._hosts = hosts_module.resolve_hosts(
                self.config, local_workdir=self.config.workdir)
        return hosts_module.host_for(instance, self._hosts)

    def _log(self, task):
        return execution.log_path(task.instance.work_dir, step_key(
            task.instance, task.phase))

    def prepare(self):
        """Create directories, resolve the instances and the skip policy."""
        self.runner.prepare()
        self._structure = self.config.structure
        instances = self.runner.selected_instances()
        self._skipped, blocking = self.runner.plan_skips(instances)
        runnable = [item for item in instances if item.id not in self._skipped]
        self._tasks, self._planned = self.runner.build_tasks(runnable)
        for instance in instances:
            if instance.id in self._skipped:
                for phase in self._planned_phases(instance):
                    self.state.mark_step(
                        step_key(instance, phase), "skipped",
                        pipeline=str(self.config.path),
                        instance=instance.id, stage=instance.stage,
                        phase=phase, host=self._host(instance).name,
                        detail={"reason": self._skipped[instance.id]})
        self._sync_state()
        return self._tasks

    def _planned_phases(self, instance):
        planned = self._planned.get(instance.id)
        if planned is not None:
            return planned
        return self.runner._planned_phases(instance)

    # ------------------------------------------------------------------
    # Config change detection
    # ------------------------------------------------------------------
    def _sync_state(self):
        """Reset the steps whose settings changed since the last run."""
        for task in self._tasks:
            host = self._host(task.instance)
            digest = config_hash(task.instance, task.phase, task.scripts,
                                 host.name, self._structure)
            self._hashes[task.key] = digest
            row = self.state.get_step(task.label)
            if row is None:
                continue
            if self.fresh or (row.get("config_hash") not in (None, digest)):
                reason = ("--fresh was requested" if self.fresh
                          else "its settings changed")
                print(f"  resetting step {task.label}: {reason}")
                self.state.reset_step(task.label, keep_record=True)
                execution.clear_manifest(task.instance.work_dir)
                self._jobs.pop(task.label, None)

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------
    def print_plan(self):
        print("\n=== Pipeline plan ===")
        for task in self._tasks:
            host = self._host(task.instance)
            phase_mode = self._mode(task)
            launcher = execution.launcher_for(
                {}, host, self._header(task.instance, task.phase))
            flags = []
            if task.wait:
                flags.append(f"calculations: {launcher}")
            if task.depends:
                flags.append("after " + ", ".join(
                    f"{i}.{p}" if p else f"{i} (data)"
                    for i, p in task.depends))
            suffix = f" [{'; '.join(flags)}]" if flags else ""
            print(f"  {task.label:28s} {phase_mode:9s} @{host.name:10s} "
                  f"{', '.join(task.scripts)}{suffix}")
            if task.instance.execution.get("max_parallel", 1) > 1:
                print(f"      max_parallel = "
                      f"{task.instance.execution['max_parallel']}")

    def report(self):
        """Human readable state summary (also used by ``--status``)."""
        lines = []
        for task in self._tasks:
            row = self.state.get_step(task.label) or {}
            done, total = self.state.items_progress(task.label)
            extra = f"  {done}/{total} calculations" if total else ""
            host = row.get("host") or self._host(task.instance).name
            lines.append(f"  {task.label:28s} {row.get('status', 'pending'):9s}"
                         f" @{host}{extra}")
            for job in self.state.tasks(
                    task.label, status=("pending", "submitted", "running")):
                lines.append(f"      {job['kind'] or '-':8s} "
                             f"id={job['external_id'] or '-'}")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------------
    def run(self):
        if not self._tasks:
            print("Nothing to do: every instance was skipped.")
            return 0
        if self._hashes is None:
            self._hashes = {}
        if self.dry_run:
            print("\n(dry-run) no command is executed; the plan only")
            self.print_plan()
            return 0
        self.print_plan()
        self.state.set_meta("config_path", str(self.config.path))
        self.state.set_meta("updated_at", time.time())
        while True:
            outcome = self._pass()
            self.state.set_meta("updated_at", time.time())
            if outcome == "done":
                self.state.close()
                print("\nPipeline finished: every step is complete.")
                return 0
            if outcome == "blocked":
                raise ConfigError(
                    "Cannot schedule the pipeline: circular dependency or a "
                    "missing prerequisite.\n" + self.report())
            if outcome == "waiting":
                if self.mode in ("wait", "daemon"):
                    pending = self._pending_steps()
                    print(f"\n  waiting for: {', '.join(pending)} "
                          f"(polling every {self.poll_interval:.0f}s)")
                    time.sleep(self.poll_interval)
                    continue
                print("\n=== Stopped: the pipeline has to wait for external "
                      "programs ===")
                print(self.report())
                print("\nRun the same command again to continue "
                      "(or use --daemon / --wait).")
                self.state.close()
                return 0
        return 0

    def _pending_steps(self):
        return [task.label for task in self._tasks
                if not self.state.is_finished(task.label)]

    # ------------------------------------------------------------------
    # One scheduling pass
    # ------------------------------------------------------------------
    def _pass(self):
        done_keys = {task.key for task in self._tasks
                     if self.state.is_finished(task.label)}
        waiting, progressed, failed = False, False, []
        for task in self._tasks:
            if task.key in done_keys:
                continue
            if not self._dependencies_met(task, done_keys):
                continue
            status = (self.state.get_step(task.label) or {}).get("status")
            if status in ("running", "waiting"):
                outcome = self._supervise(task)
            else:
                outcome = self._start(task)
            if outcome == "done":
                done_keys.add(task.key)
                progressed = True
            elif outcome == "waiting":
                waiting = True
            elif outcome == "failed":
                failed.append(task.label)
                progressed = True
            else:
                progressed = True
        if all(self.state.is_finished(task.label) for task in self._tasks):
            return "done"
        if failed:
            raise ConfigError(
                "These steps failed: " + ", ".join(failed)
                + "\nSee the logs under " + str(
                    Path(self.config.workdir) / execution.LOG_DIR))
        if progressed:
            return "progress"
        if waiting:
            return "waiting"
        return "blocked"

    def _dependencies_met(self, task, done_keys):
        for dep_id, dep_phase in task.depends:
            if dep_phase is not None and (dep_id, dep_phase) in {
                    other.key for other in self._tasks}:
                if (dep_id, dep_phase) not in done_keys:
                    return False
            else:
                instance = next((item for item in self.config.instances
                                 if item.id == dep_id), None)
                if instance is None or not self.runner._stage_done(instance):
                    return False
        return True

    # ------------------------------------------------------------------
    # Starting a step
    # ------------------------------------------------------------------
    def _mode(self, task, phase=None):
        phase = phase or task.phase
        host = self._host(task.instance)
        return execution.resolve_mode(
            task.instance, phase, host=host,
            override=self.execution_override)

    def _start(self, task):
        instance = task.instance
        host = self._host(instance)
        mode = self._mode(task)
        print(f"\n=== {task.label} ({instance.stage}, {instance.calculator}) "
              f"===")
        print(f"  dir:   {instance.work_dir}")
        print(f"  host:  {host.describe()}")
        print(f"  mode:  {mode}")

        if task.wait_before and not self.runner._wait_before(
                task, run_selected="run" in self._planned.get(instance.id, ())):
            print(f"  skipping {task.label}")
            self.state.mark_step(task.label, "skipped",
                                 detail={"reason": "no Monte Carlo output"})
            return "done"

        if instance.id not in self._input_written:
            input_file = self.runner.write_stage_input(instance)
            print(f"  input: {input_file}")
            self._input_written.add(instance.id)

        self.state.mark_step(
            task.label, "running", pipeline=str(self.config.path),
            instance=instance.id, stage=instance.stage, phase=task.phase,
            host=host.name, config_hash=self._hashes.get(task.key),
            detail={"mode": mode, "scripts": list(task.scripts),
                    "max_parallel": instance.max_parallel})
        if task.wait:
            return self._start_submit(task)
        return self._start_sync(task)

    def _step_env(self, instance, submit=False):
        env = self.runner.environment()
        if submit:
            env["AUTOMAG_AUTO_SUBMIT"] = "0"
            env[execution.MANIFEST_ENV] = str(
                execution.manifest_path(instance.work_dir))
        return env

    # -- synchronous phases ---------------------------------------------
    def _start_sync(self, task):
        instance = task.instance
        host = self._host(instance)
        mode = self._mode(task)
        if mode == "scheduler":
            return self._submit_step_job(task)

        code = 0
        for script in task.scripts:
            code = self._run_script(task, script, mode)
            if code != 0:
                break
        if code != 0:
            self.state.mark_step(task.label, "failed",
                                 detail={"returncode": code})
            print(f"  ERROR: {task.label} failed with code {code} "
                  f"(log: {self._log(task)})")
            return "failed"
        self.state.mark_step(task.label, "done", detail={"returncode": 0})
        return "done"

    def _run_script(self, task, script, mode):
        instance = task.instance
        host = self._host(instance)
        directory = STAGES[instance.stage].directory
        relative = f"{directory}/{script}"
        script_path = code_root() / relative
        if not script_path.is_file():
            raise FileNotFoundError(f"Stage script not found: {script_path}")
        print(f"  running: {relative} ({mode})")
        env = self._step_env(instance, submit=task.wait)
        if mode == "inprocess":
            from common import SubmitManual as submit_manual

            reset = getattr(submit_manual, "reset_submission_state", None)
            if reset is not None:
                reset()
            return execution.run_inprocess(script_path, instance.work_dir, env)
        command = execution.automag_command(host, relative)
        step_env = {key: value for key, value in env.items()
                    if key.startswith("AUTOMAG_") or key == "PYTHONPATH"}
        code, _, _ = execution.run_sync(host, command, instance.work_dir,
                                        log=self._log(task), env=step_env)
        return code

    # -- submit phases ---------------------------------------------------
    def _start_submit(self, task):
        instance = task.instance
        work_dir = instance.work_dir
        # The submit script regenerates the whole manifest: drop the previous
        # one so that stale or duplicated entries are never submitted.
        execution.clear_manifest(work_dir)
        self._jobs.pop(task.label, None)
        code = 0
        for script in task.scripts:
            code = self._run_script(task, script, self._mode(task))
            if code != 0:
                break
        if code != 0:
            self.state.mark_step(task.label, "failed",
                                 detail={"returncode": code})
            print(f"  ERROR: {task.label} failed with code {code} "
                  f"(log: {self._log(task)})")
            return "failed"
        jobs = self._step_jobs(task)
        jobs.register()
        self.state.mark_step(task.label, "waiting",
                             detail={"mode": self._mode(task),
                                     "calculations": jobs.count(),
                                     "max_parallel": instance.max_parallel})
        return self._supervise_submit(task)

    def _step_jobs(self, task):
        key = task.label
        jobs = self._jobs.get(key)
        if jobs is None:
            instance = task.instance
            row = self.state.get_step(key) or {}
            jobs = execution.StepJobs(
                self.state, key, self._host(instance), instance.work_dir,
                header=self._header(instance, task.phase),
                calculator=instance.calculator,
                max_parallel=instance.max_parallel,
                max_attempts=instance.max_attempts,
                dry_run=self.dry_run,
                min_mtime=row.get("reset_at"),
                log=self._log(task))
            jobs.register()
            self._jobs[key] = jobs
        return jobs

    @staticmethod
    def _header(instance, phase):
        plan = instance.phase_execution(phase)
        return plan.get("header") or instance.jobheader

    def _jobs_row(self, task):
        return self.state.get_step(task.label)

    def _supervise_submit(self, task):
        jobs = self._step_jobs(task)
        if jobs.count() == 0:
            return self._supervise_fallback(task)
        if self.mode in ("wait", "daemon"):
            finished = jobs.poll()
            if finished:
                return self._finish_submit(task, jobs)
            return "waiting"
        # one-shot: submit what fits into max_parallel, record it and leave
        jobs.pump()
        if all(item["status"] == "done"
               for item in self.state.items(task.label)):
            return self._finish_submit(task, jobs)
        return "waiting"

    def _finish_submit(self, task, jobs):
        done, total = jobs.progress()
        missing = [item for item in self.state.items(task.label)
                   if item["status"] != "done"]
        detail = {"calculations": total, "done": done}
        if missing:
            print(f"  WARNING: {len(missing)} calculation(s) of {task.label} "
                  "have no result; continuing with the available ones")
            detail["missing"] = len(missing)
        self.state.mark_step(task.label, "done", detail=detail)
        return "done"

    def _supervise_fallback(self, task):
        """No manifest: use the result files (old ``--wait`` behaviour)."""
        instance = task.instance
        stage = task.stage
        phase = task.phase
        if not task.wait:
            return "done"
        finished, expected, label = readiness.task_progress(
            stage, phase, instance.config, self.config.workdir,
            self.runner.formula(), stage_work_dir=instance.work_dir)
        if expected is None:
            print(f"  WARNING: no task manifest and no way to count the "
                  f"results of {task.label}; treating it as complete")
            self.state.mark_step(task.label, "done",
                                 detail={"source": "no manifest"})
            return "done"
        if finished >= expected:
            self.state.mark_step(task.label, "done",
                                 detail={"source": "result files"})
            return "done"
        if self.mode in ("wait", "daemon"):
            readiness.wait_for_task(
                stage, phase, instance.config, self.config.workdir,
                self.runner.formula(), timeout=self.runner.wait_timeout,
                interval=self.runner.wait_interval,
                stable_seconds=self.runner.wait_stable_seconds,
                stage_work_dir=instance.work_dir)
            self.state.mark_step(task.label, "done",
                                 detail={"source": "result files"})
            return "done"
        print(f"  {label}: {finished}/{expected} results so far")
        return "waiting"

    # -- step scripts as batch jobs --------------------------------------
    def _submit_step_job(self, task):
        instance = task.instance
        host = self._host(instance)
        relative = (f"{STAGES[instance.stage].directory}/{task.scripts[0]}"
                    if len(task.scripts) == 1 else None)
        commands = []
        for script in task.scripts:
            commands.append(execution.automag_command(
                host, f"{STAGES[instance.stage].directory}/{script}"))
        command = " && ".join(commands)
        header = self._header(instance, task.phase)
        script, sentinel = execution.step_wrapper_script(
            host, instance.work_dir, task.label, command, header=header)
        launcher = ("sbatch" if host.scheduler == "slurm" or
                    execution.scheduler_directives(header) == "slurm"
                    else "qsub" if host.scheduler == "pbs" or
                    execution.scheduler_directives(header) == "pbs"
                    else "sh")
        task_key = self.state.add_task(
            task.label, host=host.name,
            kind="slurm" if launcher == "sbatch" else
                 "pbs" if launcher == "qsub" else "process",
            workdir=str(instance.work_dir),
            command=f"{launcher} {script}", items=["__step__"],
            launcher=launcher)
        try:
            if launcher == "sh":
                external_id = host.start_background(
                    f"sh {shlex.quote(str(script))}",
                    cwd=instance.work_dir, log_path=str(sentinel))
            else:
                external_id, _ = host.submit(script, cwd=instance.work_dir,
                                             launcher=launcher)
        except (execution.HostError, OSError) as error:
            self.state.update_task(task_key, status="failed",
                                   detail={"error": str(error)})
            self.state.mark_step(task.label, "failed",
                                 detail={"error": str(error)})
            return "failed"
        self.state.update_task(task_key, status="running",
                               external_id=external_id)
        self.state.mark_step(task.label, "waiting",
                             detail={"mode": "scheduler",
                                     "batch": str(script)})
        return "waiting"

    def _supervise_step_job(self, task):
        instance = task.instance
        host = self._host(instance)
        sentinel = execution.step_sentinel(instance.work_dir, task.label)
        if sentinel.exists() if host.is_local else host.exists(sentinel):
            for job in self.state.tasks(task.label,
                                        status=("pending", "submitted", "running")):
                self.state.update_task(job["task_key"], status="finished")
            self.state.mark_step(task.label, "done",
                                 detail={"mode": "scheduler"})
            print(f"  {task.label}: batch job finished")
            return "done"
        for job in self.state.tasks(task.label,
                                    status=("pending", "submitted", "running")):
            state = host.job_state(job["external_id"], kind=job["kind"])
            if state["state"] == "failed":
                self.state.update_task(job["task_key"], status="failed")
                self.state.mark_step(task.label, "failed",
                                     detail={"scheduler": state["detail"]})
                print(f"  ERROR: batch job of {task.label} failed")
                return "failed"
            if state["state"] in ("finished", "unknown"):
                # The scheduler lost the job without the sentinel: treat it
                # as a failure, the step will be retried on the next run.
                self.state.update_task(job["task_key"], status="finished")
                self.state.mark_step(
                    task.label, "failed",
                    detail={"reason": "the batch job ended without writing "
                                      f"{sentinel}"})
                print(f"  ERROR: batch job of {task.label} ended without "
                      f"producing {sentinel.name}")
                return "failed"
        return "waiting"

    # ------------------------------------------------------------------
    # Supervision of a step that is already running
    # ------------------------------------------------------------------
    def _supervise(self, task):
        if self._timeout_reached(task):
            raise TimeoutError(
                f"Timed out waiting for {task.label} "
                f"({self.runner.wait_timeout:.0f}s).")
        status = (self.state.get_step(task.label) or {}).get("status")
        if task.wait:
            if status == "running":
                # automag stopped while the input files were being written
                return self._start_submit(task)
            return self._supervise_submit(task)
        open_jobs = self.state.tasks(
            task.label, status=("pending", "submitted", "running"))
        if status == "running" and open_jobs:
            return self._supervise_step_job(task)
        if status == "running":
            # the step script was interrupted: run it again
            return self._start(task)
        return self._supervise_step_job(task)

    def _timeout_reached(self, task):
        row = self.state.get_step(task.label) or {}
        started = row.get("started_at")
        if not started or not self.runner.wait_timeout:
            return False
        return time.time() - float(started) >= self.runner.wait_timeout
