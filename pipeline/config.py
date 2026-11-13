"""Config loading for the automag pipeline.

Supports legacy ``run`` lists (one instance per stage) and instance lists
(``run``/``stages`` entries as dicts) with per-instance calculators.
All paths resolve relative to the launch directory (``Path.cwd()``).

The controlling file may be ``input_automag.py``, ``input_automag.yaml`` or
``input_automag.json``; they are tried in that order (see
:mod:`common.input_loader`).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from common.input_loader import load_input, resolve_input
from common.paths import code_root
from pipeline.stages import STAGES, STAGE_ORDER

#: Stem of the controlling file searched for in the launch directory.
CONFIG_BASE = "input_automag"

COMMON_KEYS = (
    "calculator", "use_fireworks", "jobheader", "calculator_command",
    "environment_activate", "environment_deactivate",
    "parallel_over_configurations", "exec_command", "auto_submit",
    "external_codes",
)

#: ``header`` is accepted as an alias of ``jobheader`` (both name the batch
#: script preamble that decides whether a step is a SLURM/PBS job).
HEADER_ALIASES = {"header": "jobheader"}

CONTROL_KEYS = {
    "workdir", "structure", "struct_suffix", "run", "stages", "submit",
    "common", "calculators", "calculator_plugins", "steps", "enabled",
    "description", "wait_timeout", "wait_interval", "wait_stable_seconds",
    "wait_for_monte_carlo_output",
    # execution / scheduling
    "hosts", "execution", "state_file", "poll_interval", "max_attempts",
    "max_parallel_tasks", "parallel", "header",
}

INSTANCE_KEYS = {
    "id", "name", "stage", "calculator", "work_dir", "depends_on",
    "on_missing", "on_unsupported", "reuse", "coll_calculator",
    "coll_instance", "description", "steps", "enabled",
    # per-instance input geometry (a file in <workdir>/geometries, e.g. the
    # relaxed structure written by 2_relax) and per-instance result suffix
    "structure", "struct_suffix",
    # execution / scheduling
    "host", "execution", "max_parallel", "max_parallel_tasks", "parallel",
    "max_attempts", "header",
}

#: Execution modes understood by :mod:`pipeline.execution`.
EXECUTION_MODES = ("auto", "inprocess", "process", "scheduler")


class ConfigError(ValueError):
    pass


@dataclass
class StageInstance:
    id: str
    stage: str
    calculator: str
    config: dict
    work_dir: Path
    coll_calculator: str | None = None
    coll_instance: str | None = None
    depends_on: tuple = ()
    on_missing: str = "error"
    on_unsupported: str = "error"
    reuse: bool = False
    description: str = ""
    #: Host name (or inline configuration) the steps of this instance run on.
    host: object = None
    #: Resolved execution plan: ``{"mode": ..., "max_parallel": ...,
    #: "max_attempts": ..., "phases": {...}}``.
    execution: dict = field(default_factory=dict)

    @property
    def label(self):
        return f"{self.id}[{self.calculator}]"

    @property
    def jobheader(self):
        """Batch-script preamble of this instance (``header`` is an alias)."""
        return self.config.get("jobheader") or self.config.get("header") or ""

    @property
    def max_parallel(self):
        return int(self.execution.get("max_parallel") or 1)

    @property
    def max_attempts(self):
        return int(self.execution.get("max_attempts") or 3)

    def phase_execution(self, phase):
        """Execution settings of one phase (instance-wide + per-phase)."""
        plan = dict(self.execution)
        phases = dict(plan.pop("phases", {}) or {})
        plan.update(phases.get(phase) or {})
        return plan


@dataclass
class PipelineConfig:
    path: Path
    base_dir: Path
    workdir: Path
    structure: Path
    struct_suffix: str = ""
    run: list = field(default_factory=list)
    common: dict = field(default_factory=dict)
    stages: dict = field(default_factory=dict)
    calculator_configs: dict = field(default_factory=dict)
    calculator_plugins: list = field(default_factory=list)
    instances: list = field(default_factory=list)
    raw: dict = field(default_factory=dict)
    wait_timeout: float = 30 * 24 * 3600
    wait_interval: float = 60.0
    wait_stable_seconds: float = 600.0
    wait_for_monte_carlo_output: bool = False
    #: Named execution hosts (``PIPELINE["hosts"]``).
    hosts: dict = field(default_factory=dict)
    #: Default execution plan shared by every stage instance.
    execution: dict = field(default_factory=dict)
    #: SQLite file holding the pipeline state (None = next to the results).
    state_file: Path | None = None
    poll_interval: float = 30.0
    max_attempts: int = 3

    @classmethod
    def load(cls, path=None):
        path = _resolve_config_path(path)
        namespace = load_input(path)
        raw = namespace.get("PIPELINE") or namespace.get("pipeline")
        if raw is None:
            raw = {k: namespace[k] for k in (*CONTROL_KEYS, *STAGE_ORDER)
                   if k in namespace}
        if not isinstance(raw, dict):
            raise ConfigError(f"{path} must define a PIPELINE mapping.")

        base_dir = path.parent
        wd_raw = raw.get("workdir")
        if wd_raw is None:
            workdir = Path.cwd().resolve()
        else:
            workdir = Path(wd_raw).expanduser()
            if not workdir.is_absolute():
                workdir = (Path.cwd() / workdir).resolve()
        _check_workdir(workdir)
        workdir.mkdir(parents=True, exist_ok=True)

        if not raw.get("structure"):
            raise ConfigError("'structure' is required in the input file.")
        structure = _resolve_structure(raw["structure"], base_dir, workdir)

        common = {}
        for key in ("submit", "common"):
            if isinstance(raw.get(key), dict):
                common.update(raw[key])
        for key in COMMON_KEYS:
            if key in raw:
                common[key] = raw[key]
        if "header" in raw and "jobheader" not in common:
            common["jobheader"] = raw["header"]
        common.setdefault("calculator", "vasp")
        _apply_header_alias(common)

        calc_cfg = raw.get("calculators", {}) or {}
        if not isinstance(calc_cfg, dict):
            raise ConfigError("'calculators' must be a dictionary.")
        calc_cfg = {str(k): dict(v or {}) for k, v in calc_cfg.items()}

        stage_defaults = {}
        stage_entries = None
        stages_raw = raw.get("stages")
        if isinstance(stages_raw, (list, tuple)):
            if all(isinstance(i, dict) for i in stages_raw):
                stage_entries = list(stages_raw)
            elif all(isinstance(i, str) for i in stages_raw):
                stage_defaults = {n: {} for n in stages_raw}
            else:
                raise ConfigError("'stages' list must contain dicts or names.")
        elif isinstance(stages_raw, dict):
            stage_defaults = {n: dict(c or {}) for n, c in stages_raw.items()}
        elif stages_raw is not None:
            raise ConfigError("'stages' must be a list or a dictionary.")
        for name in STAGE_ORDER:
            if isinstance(raw.get(name), dict):
                stage_defaults.setdefault(name, {}).update(raw[name])

        run_raw = raw.get("run")
        entries = None
        if isinstance(run_raw, (list, tuple)) and any(
                isinstance(i, dict) for i in run_raw):
            entries = list(run_raw)
        elif stage_entries is not None:
            entries = stage_entries

        legacy = entries is None
        if legacy:
            if run_raw is None:
                run = [n for n in STAGE_ORDER
                       if n in stage_defaults
                       and stage_defaults[n].get("enabled", True)]
            elif run_raw == "all":
                run = list(STAGE_ORDER)
            elif isinstance(run_raw, str):
                run = [i.strip() for i in run_raw.split(",") if i.strip()]
            elif isinstance(run_raw, (list, tuple)):
                run = [str(i) for i in run_raw]
            else:
                run = []
            unknown = [n for n in run if n not in STAGE_ORDER]
            if unknown:
                raise ConfigError(f"Unknown stage(s) {unknown}.")
            entries = [{"stage": n, "id": n,
                        "calculator": common.get("calculator", "vasp"),
                        "_legacy": True} for n in run]
        else:
            run = [str(i.get("stage") or i.get("name"))
                   for i in entries if isinstance(i, dict)]

        execution_default = normalise_execution(
            raw.get("execution") or raw.get("parallel") or {})
        max_parallel_default = int(
            raw.get("max_parallel_tasks")
            or execution_default.get("max_parallel")
            or 1)
        max_attempts_default = int(
            raw.get("max_attempts") or execution_default.get("max_attempts") or 3)

        instances = _build_instances(
            entries, common, calc_cfg, stage_defaults, workdir, structure,
            raw.get("struct_suffix", ""), legacy,
            execution_default=execution_default,
            max_parallel_default=max_parallel_default,
            max_attempts_default=max_attempts_default)
        _resolve_dependencies(instances)

        plugins = raw.get("calculator_plugins") or []
        if isinstance(plugins, str):
            plugins = [i.strip() for i in plugins.split(",") if i.strip()]

        state_file = raw.get("state_file")
        if state_file is not None:
            state_file = Path(state_file).expanduser()
            if not state_file.is_absolute():
                state_file = (workdir / state_file).resolve()

        return cls(
            path=path, base_dir=base_dir, workdir=workdir,
            structure=structure, struct_suffix=raw.get("struct_suffix", ""),
            run=run, common=common, stages=stage_defaults,
            calculator_configs=calc_cfg, calculator_plugins=list(plugins),
            instances=instances, raw=raw,
            wait_timeout=float(raw.get("wait_timeout", 30 * 24 * 3600)),
            wait_interval=float(raw.get("wait_interval", 60.0)),
            wait_stable_seconds=float(raw.get("wait_stable_seconds", 600.0)),
            wait_for_monte_carlo_output=bool(
                raw.get("wait_for_monte_carlo_output", False)),
            hosts=dict(raw.get("hosts") or {}),
            execution=execution_default,
            state_file=state_file,
            poll_interval=float(raw.get(
                "poll_interval", execution_default.get("poll_interval", 30.0))),
            max_attempts=max_attempts_default,
        )


def _resolve_config_path(path=None):
    """Locate the controlling file (``input_automag.py|yaml|json``).

    With no input file at all, a directory that holds nothing but a structure
    gets a default ``input_automag.py`` (see :mod:`pipeline.bootstrap`).
    """
    try:
        if path is None:
            try:
                return resolve_input(base=CONFIG_BASE)
            except FileNotFoundError:
                from pipeline.bootstrap import ensure_config

                if ensure_config(Path.cwd()) is None:
                    raise
                return resolve_input(base=CONFIG_BASE)
        candidate = Path(path).expanduser()
        if not candidate.is_absolute():
            candidate = Path.cwd() / candidate
        if candidate.is_dir():
            return resolve_input(base=CONFIG_BASE,
                                 path=candidate / CONFIG_BASE)
        return resolve_input(base=CONFIG_BASE, path=candidate)
    except FileNotFoundError as error:
        raise ConfigError(str(error)) from error


def _build_instances(entries, common, calc_cfg, stage_defaults, workdir,
                     structure, suffix, legacy, execution_default=None,
                     max_parallel_default=1, max_attempts_default=3):
    counts = Counter()
    for entry in entries:
        if isinstance(entry, dict) and entry.get("enabled", True):
            name = entry.get("stage") or entry.get("name")
            if name in STAGE_ORDER:
                counts[name] += 1

    instances, seen = [], set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ConfigError(f"Stage instance must be a dict: {entry!r}")
        stage = entry.get("stage") or entry.get("name")
        if stage not in STAGE_ORDER:
            raise ConfigError(f"Unknown stage '{stage}'.")
        if not entry.get("enabled", True):
            continue
        calc = entry.get("calculator") or common.get("calculator", "vasp")
        iid = entry.get("id") or entry.get("name") or stage
        if iid in seen:
            raise ConfigError(f"Duplicate stage instance id '{iid}'.")
        seen.add(iid)

        cfg = dict(common)
        cfg.update(calc_cfg.get(calc, {}))
        cfg.update(stage_defaults.get(stage, {}))
        for key, value in entry.items():
            if key not in INSTANCE_KEYS and key != "_legacy":
                cfg[key] = value
        if "header" in entry and "jobheader" not in entry:
            # 'header' is a per-instance alias of the batch preamble
            cfg["jobheader"] = entry["header"]
        _apply_header_alias(cfg)
        cfg["calculator"] = calc
        # A second pass over a relaxed geometry (or over a different input
        # structure) overrides the global 'structure'.  The file is looked up
        # in <workdir>/geometries at run time, because it may be produced by an
        # earlier stage (2_relax writes <name>_relaxed.vasp there).
        own_structure = entry.get("structure")
        cfg["poscar_file"] = (Path(own_structure).name if own_structure
                              else structure.name)
        # Distinguishes the results of two instances of the same stage with the
        # same calculator, e.g. CalcFold/Fe4O6 and CalcFold/Fe4O6_relaxed.
        cfg["struct_suffix"] = entry.get("struct_suffix", suffix)
        if "steps" in entry:
            cfg["steps"] = entry["steps"]

        wd_raw = entry.get("work_dir")
        if wd_raw:
            wd = Path(wd_raw).expanduser()
            if not wd.is_absolute():
                wd = workdir / wd
        else:
            base = workdir / (STAGES[stage].work_dir
                              or STAGES[stage].directory)
            wd = base if counts[stage] == 1 else base / iid

        coll_calc = entry.get("coll_calculator")
        if coll_calc is None and stage in ("relax", "monte_carlo", "mae"):
            coll_calc = calc
        cfg["coll_calculator"] = coll_calc or calc
        deps = entry.get("depends_on") or ()
        if isinstance(deps, str):
            deps = (deps,)

        stage_cfg = stage_defaults.get(stage, {}) or {}
        plan = normalise_execution(execution_default)
        _merge_execution(plan, stage_cfg.get("execution"))
        _merge_execution(plan, entry.get("execution"))
        parallel = (entry.get("max_parallel")
                    or entry.get("max_parallel_tasks")
                    or entry.get("parallel")
                    or stage_cfg.get("max_parallel"))
        if parallel is not None:
            plan["max_parallel"] = max(1, int(parallel))
        elif not plan.get("max_parallel"):
            plan["max_parallel"] = max(1, int(max_parallel_default))
        attempts = entry.get("max_attempts") or stage_cfg.get("max_attempts")
        plan["max_attempts"] = max(
            1, int(attempts or plan.get("max_attempts") or max_attempts_default))
        plan.setdefault("mode", "auto")

        instances.append(StageInstance(
            id=iid, stage=stage, calculator=calc, config=cfg, work_dir=wd,
            coll_calculator=coll_calc,
            coll_instance=entry.get("coll_instance"),
            depends_on=tuple(deps),
            on_missing=entry.get("on_missing", "error"),
            on_unsupported=entry.get(
                "on_unsupported", "error" if legacy else "skip"),
            reuse=bool(entry.get(
                "reuse", (not legacy) and stage in ("monte_carlo", "mae"))),
            description=entry.get("description") or STAGES[stage].description,
            host=entry.get("host") or stage_cfg.get("host"),
            execution=plan,
        ))
    return instances


def _apply_header_alias(settings):
    """Make ``header`` available as ``jobheader`` (the stage scripts' name)."""
    if not isinstance(settings, dict):
        return settings
    for alias, canonical in HEADER_ALIASES.items():
        if alias in settings and canonical not in settings:
            settings[canonical] = settings[alias]
    return settings


def normalise_execution(value):
    """Turn an ``execution`` setting into a dictionary.

    Accepted forms::

        "process"                      # just the mode
        {"mode": "scheduler"}          # explicit dictionary
        {"max_parallel": 8}
        {"phases": {"submit": {"mode": "inprocess"}}}

    The returned dictionary always has ``mode`` (default ``"auto"``).
    """
    if value is None or value == {}:
        return {"mode": "auto"}
    if isinstance(value, str):
        plan = {"mode": value.strip().lower()}
    elif isinstance(value, dict):
        plan = dict(value)
        mode = plan.get("mode")
        if isinstance(mode, str):
            plan["mode"] = mode.strip().lower()
    else:
        raise ConfigError(
            "'execution' must be a string or a dictionary, got "
            f"{type(value).__name__}.")
    mode = plan.get("mode", "auto")
    if mode not in EXECUTION_MODES:
        raise ConfigError(
            f"Unknown execution mode '{mode}'. Use one of "
            f"{', '.join(EXECUTION_MODES)}.")
    plan["mode"] = mode
    if plan.get("phases") and not isinstance(plan["phases"], dict):
        raise ConfigError("'execution.phases' must be a dictionary.")
    return plan


def _merge_execution(plan, value):
    """Merge ``value`` into ``plan`` in place (phases merged one level deep)."""
    if value is None:
        return plan
    extra = normalise_execution(value)
    phases = extra.pop("phases", None)
    plan.update(extra)
    if phases:
        merged = dict(plan.get("phases") or {})
        for phase, settings in phases.items():
            current = dict(merged.get(phase) or {})
            current.update(normalise_execution(settings))
            merged[phase] = current
        plan["phases"] = merged
    return plan


def _resolve_dependencies(instances):
    by_id = {i.id: i for i in instances}
    for inst in instances:
        if inst.depends_on:
            missing = [d for d in inst.depends_on if d not in by_id]
            if missing:
                raise ConfigError(
                    f"'{inst.id}' depends on unknown instance(s): {missing}.")
            continue
        if inst.stage not in ("relax", "monte_carlo", "mae"):
            continue
        if inst.coll_instance:
            if inst.coll_instance not in by_id:
                raise ConfigError(
                    f"'{inst.id}' references unknown coll_instance "
                    f"'{inst.coll_instance}'.")
            inst.depends_on = (inst.coll_instance,)
            continue
        candidates = [o for o in instances
                      if o.stage == "coll"
                      and o.calculator == inst.coll_calculator]
        if len(candidates) == 1:
            inst.depends_on = (candidates[0].id,)
        elif len(candidates) > 1:
            raise ConfigError(
                f"'{inst.id}' matches several coll instances for "
                f"'{inst.coll_calculator}'; set 'coll_instance'.")


def _check_workdir(workdir):
    code = code_root().resolve()
    resolved = workdir.resolve()
    if resolved == code:
        raise ConfigError(
            "workdir resolves to the automag checkout. Run the pipeline "
            "from its own directory (or use examples/).")
    try:
        relative = resolved.relative_to(code)
    except ValueError:
        return
    if not relative.parts or relative.parts[0] != "examples":
        raise ConfigError(
            "A workdir inside the automag checkout must be under examples/.")


def _resolve_structure(value, base_dir, workdir=None):
    path = Path(value).expanduser()
    if path.is_absolute():
        candidates = [path]
    else:
        candidates = [base_dir / path, Path.cwd() / path]
        if workdir is not None:
            candidates.append(Path(workdir) / "geometries" / path)
        candidates += [
            code_root() / "geometries" / path,
            code_root() / "geometries" / path.name,
        ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise ConfigError(
        "Structure file not found. Tried: "
        + ", ".join(str(c) for c in candidates))