# automag — dense agent reference

Read §0 + routing, then only the section you need. Code is the source of truth.
History of the last big refactor (decisions, bugs already fixed, verified/not verified
state, git status): `docs/SESSION_HANDOFF.md` — consult it before reworking anything.

## 0. TL;DR

```bash
cd <dir with input_automag.py>                     # all results land here
python /path/to/automag/run_automag.py --wait      # full pipeline
python /path/to/automag/run_automag.py --dry-run --wait    # plan only
python /path/to/automag/run_automag.py --list-stages
python /path/to/automag/run_automag.py --stages coll,relax
python /path/to/automag/run_automag.py --instances coll_prophet,mc_vasp
./build.sh                                         # -> dist/run_automag (ninja)
```

- Config: `PIPELINE` dict in `input_automag.py`; may live anywhere; `workdir`
  defaults to the launch directory (`Path.cwd()`), so results land next to the
  config and the automag checkout is never written to.
- Results: `./{0_conv,1_lin_response,2_coll,2_relax,3_monte_carlo,4_mae}/`,
  `./CalcFold`, `./geometries`.
- Stages may appear **several times** with different calculators (instances),
  see §3.
- Examples: `examples/{vasp,qe,fhi-aims,pymatgen_enumeration,prophet,mixed_*,remote_*}/`
  (copy-and-run, each with six `run_pipeline*.sh`: `--wait`, `--once`,
  `--once-loop`, compiled and SLURM driver variants).  `prophet/` uses the
  Prophet-V1-Mag checkpoint instead of DFT (`PROPHET_MODEL`, `PROPHET_DEVICE`);
  `remote_vasp|remote_fhi_aims|remote_qe` run `2_coll` on `hosts["cluster"]`
  (ssh/shell transport) - `tools/remote_examples_test.sh` exercises that path.
- Do not run real DFT unless asked. Validate with §10.

### Routing: task → files to read/modify

| Task | Read/modify |
| --- | --- |
| Run/inspect pipeline | `run_automag.py`, `input_automag.py`, `pipeline/*` |
| Add/change calculator | `calculators/base.py`, `calculators/<name>.py`, `calculators/__init__.py` |
| Add/change stage or phase | `pipeline/stages.py`, `pipeline/readiness.py` |
| Change config keys/merge | `pipeline/config.py` |
| Change job scripts/inputs | `common/SubmitManual.py`, `common/write_output*.py` |
| Post-processing | `<stage>/*.py`, `common/results.py` |
| Examples | `examples/` |

Do not read `CalcFold/`, `runs/`, `geometries/`, archives unless debugging.

## 1. Repo map

| Path | Role |
| --- | --- |
| `run_automag.py` | Entry point; bootstraps `sys.path`/`AUTOMAG_CODE_PATH`; `--run-script`/`--run-helper` |
| `input_automag.py` | Documented example config |
| `pipeline/config.py` | `PipelineConfig.load`, `StageInstance`, `ConfigError`, key lists |
| `pipeline/stages.py` | `Stage` dataclass, `STAGES`, `STAGE_ORDER` (dependency graph) |
| `pipeline/runner.py` | `PipelineRunner`, `Task`, skip/reuse policies, `stage_summary` |
| `pipeline/readiness.py` | Polling predicates for `--wait` |
| `pipeline/cli.py` | CLI (`--wait`, `--stages`, `--instances`, `--steps`, `--dry-run`) |
| `calculators/` | Interface registry + VASP/QE/FHI-aims/prophet/generic |
| `common/paths.py` | `code_root()`, `data_root()`, `calcfold()`, `geometries_dir()`, `frozen()`, `helper_command()` |
| `common/magnetic.py` | `resolve_spin_values`, `resolve_magnetic_elements`, `pymatgen_configurations` (DEFAULT_MAGMOMS + both enumeration backends) |
| `common/results.py` | `2_coll` artefacts: `coll_bundle`, `available_coll_calculators`, `most_stable_state`, `find_trials_dir`, `find_state`, `load_state`, `resolve_configuration` |
| `common/input_loader.py` | `resolve_input`, `load_input`, `apply_stage_input` (`.py` → `.yaml` → `.json`) |
| `common/SubmitManual.py` | Writes calculator inputs + job scripts |
| `common/SubmitFirework.py` | Legacy FireWorks path (unmaintained; prefer manual) |
| `common/write_output.py` | VASP result writer → `CalcFold/*.txt` |
| `common/write_output_ase.py` | Generic result writer for non-VASP |
| `common/write_charges.py`, `write_magmoms_input.py`, `write_output_vasp.py`, `get_magmoms_vasp.py` | Legacy VASP helpers |
| `common/run_prophet.py` | Runs `prophet_spin.SpinCalculator` in a job dir |
| `0_conv_tests/`, `1_lin_response/`, `2_coll/`, `2_relax/`, `3_monte_carlo/`, `4_mae/` | Stage scripts |
| `examples/` | Copyable VASP/QE/FHI-aims/pymatgen/mixed examples + launchers |
| `build.ninja`, `build.sh`, `tools/` | Single-file build (PyInstaller) |

### Bootstrap and the bundled models

- `pipeline/bootstrap.py`: with no `input_automag.*` in the launch dir but a
  structure file present, writes a default `input_automag.py` (plan:
  `coll` with `prophet_spin`; **not** `4_mae` — the released model is isotropic
  and MAE would be 0) and continues.  Never overwrites an existing file.
- `models/` holds both Prophet checkpoints; `common.paths.models_dir()` /
  `model_file()` resolve them for the source tree and for the frozen bundle.
  Examples leave `model_path` unset (falls back to the bundle) and only allow
  `PROPHET_MODEL` / `PROPHET_HEAD` overrides.

## 2. Execution model

```
input_automag.py
  └─ PipelineConfig ──> instances[StageInstance] ──> PipelineRunner
        ├─ non-wait: instances in order → selected phases (scripts in order)
        └─ --wait: DAG over (instance, phase) → run task → poll → next tasks
              └─ workdir (= launch dir by default)
                   ├─ geometries/<structure basename>
                   ├─ CalcFold/<formula><suffix>/<calculator>/<state>/{singlepoint,recalc}
                   └─ <stage work_dir>[/<instance id>]/input.py + stage outputs
```

One **instance** = one stage run with one calculator. `2_coll` may therefore
appear several times in a pipeline; each instance gets its own directory:

| situation | directory |
| --- | --- |
| a single instance of that stage | `<workdir>/<Stage.work_dir>` (e.g. `2_coll`) |
| several instances of that stage | `<workdir>/<Stage.work_dir>/<instance id>` (e.g. `2_coll/coll_prophet`) |
| explicit `work_dir` in the instance | `<workdir>/<work_dir>` (absolute paths allowed) |

`Stage.work_dir` names: `0_conv`, `1_lin_response`, `2_coll`, `2_relax`,
`3_monte_carlo`, `4_mae`. `CalcFold` and `geometries` are always global to
the workdir.

Env vars:

| Var | Meaning |
| --- | --- |
| `AUTOMAG_CODE_PATH` | Repo root (set by `run_automag.py`/runner) |
| `AUTOMAG_WORKDIR` | Project data root (runner sets) |
| `AUTOMAG_PATH` | Legacy: data root fallback; runner sets `= workdir` |

`common.paths.code_root()` prefers `AUTOMAG_CODE_PATH`; `data_root()` prefers
`AUTOMAG_WORKDIR` then `AUTOMAG_PATH`. Job scripts run with cwd under
`workdir/CalcFold/...`; helper writers derive the data root from env or by
walking up to a `CalcFold` parent.

`--wait` phases: `submit` → poll (`wait_after`) → `postprocess`; `monte_carlo.run`
is synchronous; `monte_carlo.postprocess` has a `wait_before` fallback for an
externally produced `output`.

### Execution layer (`pipeline/hosts.py`, `execution.py`, `state.py`, `scheduler.py`)

The runner does not execute phases itself any more: `PipelineRunner.run_scheduled`
hands the plan to `pipeline.scheduler.Scheduler`, which is the state machine
(`once` / `wait` / `daemon`).

| module | responsibility |
| --- | --- |
| `pipeline/hosts.py` | `Host` (transport `local`/`ssh`/`shell`, scheduler `slurm`/`pbs`/`none`, `automag_home`/`automag_exe`, `python`, `environment`, `wrapper`), `resolve_hosts`, `submit`/`submit_command` (`sbatch`/`qsub`/`sh`), `classify_state` (SLURM/PBS status words → pending/running/finished/failed/unknown), `wrap`, `run_sync` |
| `pipeline/execution.py` | step scripts (`Steps`), calculation jobs (`StepJobs`), `resolve_mode`, `launcher_for`, manifest reading, marker/`outcar_finished` checks, `run_sync`, `automag_command` |
| `pipeline/state.py` | SQLite `StateStore`: `steps`, `items`, `tasks`, `meta` tables + `MIGRATIONS`; `reset_step`, `add_items`, `report` |
| `pipeline/scheduler.py` | one pass = pick steps whose dependencies are done, start/supervise/finish them; `config_hash` per step detects changed settings |

Key behaviours (all covered by `tools/execution_test.py`):

* **Step scripts**: `execution: auto` → `process` (separate `run_automag
  --run-script` process), or `inprocess` when `can_spawn_process()` says the
  machine cannot afford another process; `scheduler` is explicit only. A step on
  a non-local host always runs as a process there.
* **Calculations**: the step writes a JSON-lines manifest
  (`AUTOMAG_TASK_MANIFEST`, written by `common/SubmitManual._record_task`);
  automag reads it back and starts each entry with `launcher_for` — an explicit
  `launcher`, else the host scheduler, else the `#SBATCH`/`#PBS` directives in
  the step `header`, else a plain process. `max_parallel` bounds the concurrent
  jobs; `dependencies_finished` honours per-entry `depends_on` (this is what
  orders `4_mae` `z` → θ/φ → curve).
* **Item tracking**: every calculation is an `items` row with a marker path
  (`.automag_done` / `outcar_finished`); a task with no marker after
  `max_attempts` is left unfinished, the step continues and warns.
* **Config change**: `Scheduler` hashes `instance.config` per phase; a different
  hash calls `StateStore.reset_step` (status `pending`, `reset_at`), drops the
  step's items/tasks, and `_spec_satisfied` ignores markers older than
  `reset_at`, so the affected calculations are recomputed while untouched steps
  keep their `done` status.
* **`--once`** runs exactly one pass, records the started jobs, prints
  `=== Stopped: the pipeline has to wait for external programs ===` and exits 0;
  `--wait`/`--daemon` sleep `poll_interval` and poll (`--status`/`--report` read
  the state without doing anything).

Parameter sweeps (`common/sweep.py`): a *list* value of one of `SWEEPABLE`
(`encut`, `kpts`, `sigma`, ...) in a stage's `params` expands into one
calculation set per combination; `4_mae/1_submit.py` writes each set into its
own `mae_*_S01` directory with `parameters.json` and records the list in
`sweep_sets.json`; `2_plot_results.py` / `3_submit.py` / `4_plot_results.py`
iterate over exactly those directories and write
`MAE_convergence_<conf>.txt` / `MAE_curve_convergence_<conf>.txt`. Without a
manifest the scripts fall back to globbing (historical, single-set behaviour,
no `_S01` suffix).

## 3. Config schema (`input_automag.py`)

**Format**: the controlling file is looked up as `input_automag.py` →
`input_automag.yaml`/`.yml` → `input_automag.json` (priority in that order,
via `common.input_loader.resolve_input`). YAML/JSON may wrap the dict in a
`PIPELINE:` key or be the dict itself; a flat mapping works too. The same
fallback applies to the per-stage `input.*` files (`input.py` → `input.yaml` →
`input.json`), so `common.input_loader.apply_stage_input(globals())` replaces
the old `from input import *`. In YAML/JSON a `range(...)` becomes a plain list
(the stage code only iterates).

Two equivalent forms. Legacy form — one instance per stage:

```python
PIPELINE = {
  "workdir": ".",                   # default: launch dir (Path.cwd())
  "structure": "structure.vasp",     # absolute / relative / name in automag/geometries
  "struct_suffix": "",
  "run": ["conv_tests", "coll"],   # list, comma-string, or "all"
  "calculator_plugins": [],          # module paths registered in calculators
  "wait_timeout": 30*24*3600,
  "wait_interval": 60,
  "wait_stable_seconds": 600,
  "wait_for_monte_carlo_output": False,
  "submit": { ... },                 # shared job settings
  "coll": { ... }, "relax": { ... }, # per-stage dicts
}
```

Instance form — any order, repeats allowed (`run` or `stages` as a list of dicts):

```python
PIPELINE = {
  "workdir": ".", "structure": "structure.vasp",
  "run": [
    {"id": "coll_prophet", "stage": "coll", "calculator": "prophet_spin",
     "params": PROPHET_PARAMS, "supercell_size": 1, "spin_values": {...}},
    {"id": "relax_vasp", "stage": "relax", "calculator": "vasp",
     "coll_calculator": "prophet_spin", "on_missing": "skip",
     "params": VASP_PARAMS},
    {"id": "mc_vasp", "stage": "monte_carlo", "calculator": "vasp",
     "coll_calculator": "vasp", "on_missing": "skip", "configuration": "lowest"},
  ],
}
```

Instance keys: `id`, `stage`, `calculator`, `work_dir`, `depends_on`,
`coll_calculator`, `coll_instance`, `on_missing`, `on_unsupported`, `reuse`,
`enabled`, `steps`, `description`. **Any other key is a stage setting**
(`params`, `supercell_size`, `jobheader`, ...).

Merge order for the generated `input.py`:
`submit`/`common` → `calculators[<calc>]` → `stages[<stage>]` → instance dict
(+ always `calculator`, `poscar_file`, `struct_suffix`, `coll_calculator`).
`submit`/`common` keys: `calculator`, `use_fireworks`, `jobheader`,
`calculator_command`, `environment_activate`, `environment_deactivate`,
`parallel_over_configurations`, `auto_submit`, `exec_command`.

Policies:

| Key | Default | Meaning |
| --- | --- | --- |
| `on_missing` | `"error"` | `"skip"`: warn + skip when prerequisite results are absent |
| `on_unsupported` | `"skip"` (dict) / `"error"` (legacy) | warn + skip when the calculator cannot run the stage |
| `reuse` | `True` for `monte_carlo`/`mae` (dict form), else `False` | skip if the output already exists (warn) |

- `coll_calculator` = whose `2_coll` bundle this instance reads; defaults to the
  instance's own `calculator`. Set it to read another calculator's results.
- Auto dependency: `relax`/`monte_carlo`/`mae` depend on the single `coll`
  instance whose calculator == `coll_calculator`; with several candidates set
  `coll_instance` (or `depends_on`).
- Capability gate (`runner.CAPABILITIES`): `lin_response`→`supports_ldau`,
  `relax`→`supports_relaxation`, `mae`→`supports_noncollinear`.
- `configuration = "lowest"` in `monte_carlo`/`mae`/`relax` resolves to the most
  stable state of `coll_calculator` (`common.results.resolve_configuration`).
- `workdir` may not be the automag root; inside the checkout it must be under
  `examples/`.
- Values must be `pprint`-serialisable (lists/ranges/dicts/str; no numpy arrays).
- Generated `input.py` preserves `# LINE ADDED BY THE SCRIPT ...` blocks.
- Flat module variables are also accepted if `PIPELINE` is absent.

## 4. Stage graph (`pipeline/stages.py` is authoritative)

| stage | dir | default | `--wait` full | `wait_after` | depends (phase ← deps) | required keys |
| --- | --- | --- | --- | --- | --- | --- |
| `conv_tests` | `0_conv_tests` | submit | submit, postprocess | submit | postprocess ← conv_tests.submit | mode, params |
| `lin_response` | `1_lin_response` | submit | submit, postprocess | submit | postprocess ← lin_response.submit | dummy_atom, dummy_position, perturbations, params |
| `coll` | `2_coll` | submit | submit, postprocess | submit | postprocess ← coll.submit | supercell_size, params |
| `relax` | `2_relax` | submit | submit, postprocess | submit | submit ← coll.postprocess; postprocess ← relax.submit | params |
| `monte_carlo` | `3_monte_carlo` | submit | submit, run, postprocess | — | submit ← coll.postprocess; run ← monte_carlo.submit; espins ← submit; postprocess ← run | configuration, cutoff_radius, control_group_size |
| `mae` | `4_mae` | submit | submit, axes, curve, curve_postprocess | submit, curve | submit ← coll.postprocess; axes ← mae.submit; curve ← mae.axes; curve_postprocess ← mae.curve; postprocess ← mae.submit | configuration, params, Nph, Nth, N_MAE |

Readiness (`pipeline/readiness.py`):
`conv_tests.submit` line count `CalcFold/<f><s>_<mode>_<calc>.txt`;
`lin_response.submit` charges lines; `coll.submit` singlepoint lines vs
`trials_*/configurations*.txt`; `relax.submit` finished `OUTCAR`;
`mae.submit` finished `PhTh_*`; `mae.curve` finished `K_*_RtMAE_*`.
File-count checks stop early on `expected` or after `wait_stable_seconds` of no
progress.

## 5. Data contracts (do not break)

`CalcFold` files:

| File | Producer | Consumer |
| --- | --- | --- |
| `<formula><suffix>_singlepoint_<calc>.txt` | `write_output.py` / `write_output_ase.py` | `2_coll/2_plot_results.py` |
| `<formula><suffix>_<mode>_<calc>.txt` (`mode`=encut/kgrid) | same | `0_conv_tests/2_plot_results.py` |
| `<formula>_charges_<calc>.txt` (no suffix!) | `write_charges.py` | `1_lin_response/2_plot_results.py` (reads with suffix) |
| `3_monte_carlo/output` (2 cols: T, M) | VAMPIRE/ESpins or `5_run_mc.py` | `3_plot_results.py` |

Collinear line format (exact token order matters; `values[-4]` is used):

```
<state>       <calcid>  singlepoint=<convergedconverged|NONCONVERGEDNONCONVERGED>  <formula>   <spacegroup>  energy=<eV>
          initial_magmoms=[...]  final_magmoms=[...]
```

`initial_magmoms.txt` is written in `.../singlepoint/`. `2_coll` output bundle
(one per instance, hence the recursive lookup):

```
2_coll[/<instance>]/<formula><suffix>_<calc>/
    trials_<formula><suffix>/{settingNNN.vasp,configurationsNNN.txt}
    <formula><suffix>_<calc>_{states,energies,setting}NNN.txt/.vasp
    <formula><suffix>_<calc>_most_stable_configuration.txt
    <formula><suffix>_<calc>_configurations_stability_ascending.txt
```

Always read it through `common/results.py` and always pass the calculator:
`coll_bundle(root, formula, suffix, calculator)`,
`most_stable_state(..., calculator)`, `find_trials_dir(..., calculator)`,
`find_state(..., calculator=...)`, `load_state(..., calculator=...)`.
A downstream stage uses `coll_calculator` (not its own `calculator`) for these
calls and for `CalcFold/<formula><suffix>_singlepoint_<coll_calculator>.txt`.

## 6. Calculator interface (`calculators/`)

`CalculatorInterface` methods:

| Method | Purpose |
| --- | --- |
| `write_input(atoms, directory, params=None)` | Write code inputs (ASE 3.29: use `write_inputfiles` when present) |
| `run_command()` | Shell command for one calculation |
| `read_energy/read_structure/read_magnetic_moments/is_converged(dir)` | Result parsing |
| `result_writer()` | `(script, args)` appended to job script; use `common/write_output_ase.py` |
| `charges_writer()/magmoms_writer()` | VASP LDA+U/recalc helpers |

Capability flags: `supports_ldau`, `supports_noncollinear`,
`supports_relaxation`, `supports_charges`. `SubmitManual` guards VASP-only modes;
`runner.CAPABILITIES` uses the flags to skip unsupported instances early.
Currently `supports_relaxation`/`supports_noncollinear` are VASP-only (relax and
MAE emit VASP tags), so `relax`/`mae` with QE/FHI-aims/prophet are skipped.

Registry: `get_calculator(name)`, `register`, `load_plugins`. Built-ins:
`vasp`, `espresso`/`qe`, `fhi-aims`/`aims`, `prophet_spin`/`prophet`, `ase`.
`calculator = "ase:<dotted.path>"` works for any ASE calculator.

Add a calculator (no pipeline edits):
1. subclass `CalculatorInterface` in `calculators/<name>.py`;
2. implement `write_input`, `run_command`, result readers, `result_writer`;
3. `calculators.register` it and add to `load_plugins` aliases;
4. add a copyable example if it is a first-class backend;
5. validate §10.

ASE 3.29 needs explicit profiles: `EspressoProfile(command, pseudo_dir)` and
`AimsProfile(command, default_species_directory)`; already handled.

## 7. `SubmitManual` modes

| mode | VASP | other calculators |
| --- | --- | --- |
| `encut`, `kgrid` | VASP single-point | generic single-point |
| `singlepoint` | singlepoint + recalc | single-point only |
| `relax` | VASP relaxation tags | generic single-point (`params` still VASP-style) |
| `perturbations` | LDA+U/U (VASP only) | guarded `NotImplementedError` |
| `mae_theta_phi`, `mae_curve` | non-collinear MAE (VASP only) | guarded |

Non-VASP: `_write_inputs` writes only `singlepoint`, writes
`initial_magmoms.txt`, and `_write_result_writer` runs
`common/write_output_ase.py` in place (no copy). The legacy stage default
`calculator_command = "mpirun vasp_std"` is ignored for non-VASP in
`SubmitManual.__init__` so the interface command is used.

Submission: `SubmitManual.submit()` writes the inputs and then launches the
accumulated `run_all_<mode>.sh` scripts (`sbatch` / `qsub` / `sh` chosen from
the job header). Disable with `"submit": {"auto_submit": False}` — the runner
passes it as `AUTOMAG_AUTO_SUBMIT=0`. The scripts are truncated on first touch
per process, so re-running a stage never resubmits stale jobs; only commands
queued since the last flush are launched. Requires
`parallel_over_configurations = True`, otherwise nothing is submitted.

## 8. Add a stage/phase

1. Add a `Stage(...)` entry to `STAGES` and `STAGE_ORDER` in
   `pipeline/stages.py` (`steps`, `default_steps`, `full_steps`, `depends`,
   `wait_after`, `required`).
2. Add readiness for phases that submit jobs in `pipeline/readiness.py`
   (`wait_for_task` dispatch).
3. Stage scripts must: read their settings via `common.input_loader`
   (`load_input` / `apply_stage_input`), use
   `common.paths.data_root()/geometries_dir()`, write results to `CalcFold`
   or the stage dir, and reach external codes through `common.external`
   (`require`/`run`/`shell_command`) — never hard-code their paths.
4. If a new top-level config key is needed, add it to `CONTROL_KEYS` and parse
   it in `PipelineConfig.load`.
5. Validate with `--dry-run --wait`; never bypass the DAG.

## 9. Invariants / pitfalls

- Do not rename stage dirs or `CalcFold` filenames; postprocessing depends on
  them. `Stage.work_dir` (e.g. `2_coll`) is a *new* name used for the workdir
  sub-directory only — the code directories keep their old names.
- Never write into the automag checkout: results go to `workdir` (default
  `Path.cwd()`). `config._check_workdir` enforces this.
- Read `2_coll` artefacts only through `common/results.py` with an explicit
  calculator; never reconstruct the path by hand (instances add a level).
- Keep `AUTOMAG_PATH` legacy semantics; use `common.paths` for new code.
- Do not reformat or "clean up" the numerical VASP code.
- `2_coll/2_plot_results.py` uses `structure.formula` (not reduced) and a
  suffix in the result template; check every consumer before changing names.
- `write_charges.py` writes `<formula>_charges_<calc>.txt` without suffix;
  `1_lin_response/2_plot_results.py` reads with suffix — only consistent for
  empty `struct_suffix`.
- QE `starting_magnetization` is per species; AFM within one species needs
  distinct sublattice labels. FHI-aims (`initial_moment`) and prophet are
  per atom.
- Never hard-code `input.py`: use `common.input_loader.load_input` /
  `apply_stage_input` so `input.yaml` and `input.json` keep working.
- `SubmitFirework` is unmaintained and stage scripts pass `struct_suffix` it
  does not accept; prefer `use_fireworks=False`.
- `monte_carlo.run` is synchronous; `wait_for_monte_carlo_output` only matters
  if `run` is skipped.
- Preserve `# LINE ADDED BY THE SCRIPT` markers; the runner relies on them.
- `--wait` polls files, not the scheduler. File-count checks proceed after
  `wait_stable_seconds` without progress; `relax` waits for a finished
  `OUTCAR` or the timeout.
- Old per-stage `input.py` workflow still works via `AUTOMAG_PATH` fallback.
- `reuse` for `monte_carlo` checks `output` (not `model.png`); for `mae` it
  checks `MAE_curve_*.txt`; for `relax` `relaxed_*.json`.
- Skipped instances are reported as `WARNING: skipping instance ...` and their
  dependents are skipped too unless the dependency's data already exists.
- A `relax`/`monte_carlo`/`mae` instance with no matching declared `coll`
  instance falls back to checking `2_coll` data on disk for `coll_calculator`,
  and is skipped/errors according to `on_missing`.
- The stage subprocesses get the instance work directory on `PYTHONPATH`: the
  scripts do `from input import *` and `python <abs script>` does not add the
  cwd to `sys.path`.

## 10. Validation

```bash
source /home/dpoletaev/miniconda3/bin/activate automag && \
python -m py_compile run_automag.py input_automag.py pipeline/*.py \
  calculators/*.py common/*.py 2_relax/*.py 3_monte_carlo/5_run_mc.py \
  3_monte_carlo/[124]_*.py 4_mae/*.py tools/*.py

# Fast end-to-end check: no DFT code, no enumlib, ~30 s.
source /home/dpoletaev/miniconda3/bin/activate automag && \
python tools/smoke_test.py

source /home/dpoletaev/miniconda3/bin/activate automag && \
python run_automag.py --list-stages

# examples must be run from their own directory (workdir = launch dir).
# Each example also ships run_pipeline*.sh: --wait (blocking), --once and
# --once-loop (manual passes), compiled variants and a self-submitting SLURM
# driver job; they only wrap `run_automag.py [--once|--wait]`.
source /home/dpoletaev/miniconda3/bin/activate automag && \
cd examples/vasp && python ../../run_automag.py --dry-run --wait
cd examples/mixed_qe_vasp && python ../../run_automag.py --dry-run --wait
cd examples/remote_vasp && python ../../run_automag.py --list-hosts
```

`tools/smoke_test.py` fabricates a ``2_coll`` bundle in a temp dir and asserts
that `reuse`, `on_missing="skip"`, cross-calculator reads and real script
execution work. Run it after touching `pipeline/`, `common/results.py` or any
stage script - it catches what `--dry-run` cannot.

```bash
# Execution layer: --once barrier, resume, parallel limit, config-change reset,
# sbatch/qsub launcher, ssh/shell host. No DFT, no real cluster.
source /home/dpoletaev/miniconda3/bin/activate automag && python tools/execution_test.py

# Whole 4_mae chain with synthetic VASP output + the parameter sweep.
source /home/dpoletaev/miniconda3/bin/activate automag && python tools/mae_chain_test.py

# The three example launchers (blocking, one pass, loop) with a fake calculator.
source /home/dpoletaev/miniconda3/bin/activate automag && bash tools/launchers_test.sh
```

Both print a `FAILED:` list and exit non-zero on failure; `--keep` leaves the
temporary projects (under `$TMPDIR`, default
`~/.hermes/cache/scratch`) for inspection.

For calculator changes also check `calculators.available()` and write inputs in
a temp dir. Do not leave `runs/`, `examples/*/<stage>/`, `examples/*/CalcFold`
or `build/`, `dist/` behind.

## 11. Single-file build (frozen)

```bash
source /home/dpoletaev/miniconda3/bin/activate automag
conda install -c conda-forge ninja && pip install pyinstaller
./build.sh --onedir                 # ninja -> dist/run_automag/
automag_exe=dist/run_automag/run_automag; $automag_exe --list-stages

# The bundle must work with no conda, no sources and no PYTHONPATH:
env -i PATH=/usr/bin:/bin $automag_exe --dry-run

AUTOMAG_EXE=$PWD/dist/run_automag examples/vasp/run_pipeline.sh --dry-run
```

Two variants exist and are built with the same script:

| build | command | size | prophet |
| --- | --- | --- | --- |
| embedded | `--onedir --name run_automag` | 5.7 GB | `_internal/models/*.pt`, nothing to set |
| external | `--onedir --name run_automag-extmodels --no-models` | 5.4 GB | models from `PROPHET_MODEL` / `PROPHET_HEAD` |
| small | `--onedir --name run_automag-noprophet --no-prophet --no-models` | 500 MB | none (VASP/QE/FHI-aims/ASE) |

`common.paths.model_file()` tries the environment first (`PROPHET_MODEL` for
the energy model, `PROPHET_HEAD` for the head), then the bundled copy - so the
external variant works from any directory as long as the variables are set, and
the embedded one needs nothing.

`models/prophet-v1-mag.pt` and `models/prophet-v1-mag-head.pt` are copied into
the bundle (found via `common.paths.models_dir()` = `<bundle>/_internal/models`),
and `prophet_spin` + torch are included by default (`--no-prophet` drops them).
Torch needs `--collect-all torch`, `--exclude-module torch._numpy` and
`--collect-submodules ase` (for `ase.io.db`), otherwise the frozen imports fail
with `NameError: name 'name'` / `cannot import name 'detect_fake_mode'`.

Mechanics: `tools/build_executable.py` runs PyInstaller with `--add-data` for
all stage dirs + `common`/`calculators`/`pipeline`/`geometries`;
`tools/write_ninja_deps.py` regenerates `build/deps.ninja` (ninja has no glob).
In a frozen process `common.paths.code_root()` is `sys._MEIPASS` and:
- the runner starts stages as `<exe> --run-script <dir>/<script>`,
- job scripts and the prophet interface use `common.paths.helper_command()`
  → `<exe> --run-helper common/<script>`,
- paths are **relative to the bundle root** because every launch of a one-file
  build unpacks to a fresh temp dir.
`run_automag.py` implements both switches.

## 12. Magnetic configurations (`2_coll` enumeration)

`spin_values` and `magnetic_atoms` are optional everywhere. `common/magnetic.py`:

- `resolve_spin_values(structure, spin_values, magnetic_atoms, overrides=None)`
  → `{symbol: [magnitude, ...]}`, filling missing entries from
  `pymatgen.analysis.magnetism.analyzer.DEFAULT_MAGMOMS` (non-tabulated
  elements are omitted = non-magnetic). Integral moments stay `int` because
  the stage scripts format them with `{:2d}` and use them as spin
  multiplicities.
- `resolve_magnetic_elements(...)` → magnetic symbols (explicit
  `magnetic_atoms` wins, else non-zero moment).
- `structure_moments(...)`, `ensure_enumlib_on_path()`, `site_spin(site)`
  (`site.spin` does **not** exist on `PeriodicSite`; the spin lives in
  `site.specie.spin`, encoded as `Fe,spin=np.float64(-5.0)`).
- `pymatgen_configurations(...)` → `(lattices, coordinates, configurations)`
  for `enumeration_mode = "pymatgen"`
  (`MagneticStructureEnumerator`; grouping is by **lattice**, so all spin
  patterns of a cell land in one `configurationsNNN.txt`).

`2_coll/1_submit.py` selects the backend: `enumeration_mode` = `enumlib`
(default, original Wyckoff/split driver) or `pymatgen`. Both fill the same
`lattices`/`coordinates`/`configurations` lists and write the same
`trials_*/settingNNN.vasp` + `configurationsNNN.txt`, so `2_plot_results.py`,
`3_monte_carlo` and `4_mae` are unaffected. pymatgen finds `enum.x` on `PATH`,
so `ensure_enumlib_on_path()` is called before enumerating.

## 13. External codes (enumlib / VAMPIRE / ESpinS)

| code | consumer | entry keys |
| --- | --- | --- |
| `enumlib` | `2_coll/1_submit.py` | `enum`, `makestr`, `python` |
| `vampire` | `3_monte_carlo/5_run_mc.py` | `command` |
| `espins` | `4_write_espins_mcin.py`, `5_run_mc.py` | `command` |

Runtime registry: `common/external.py` (`registry`, `tool`, `available`,
`command`, `shell_command`, `require`, `run`, `summary`). Resolution order:
built-in defaults (`enum.x`/`makeStr.py`, `vampire`, `mc.x` on `PATH`) →
`code_root()/external/codes.json` (or `<exe dir>/external/codes.json` when
frozen) → `AUTOMAG_EXTERNAL_CODES` (file or dir) → `external_codes` in the
pipeline input. Entries may carry `pre_command`/`post_command` (run in one
bash, so `module`/`source` work) and `available: false`.

Build/register: `tools/build_external.py` (or `./build.sh --external`, ninja
target `external`). `--<code>-path` records an existing installation and skips
cloning/compiling; `--config FILE` does the same from JSON/YAML. A failed
compile records `available: false`, warns, and tells the user to compile
manually and pass the path — the rest of the pipeline keeps working. Recipes:
enumlib = clone `--recursive` + `make F90=gfortran` in `symlib/src` then
`enum.x makestr.x`; VAMPIRE = `make serial` → `vampire-serial`; ESpinS =
`make` → `mc.x` (needs Intel Fortran + MKL). Each build step is retried once
with the system toolchain (`system_env()`), because a conda `g++` often fails
to link against the system glibc.

## 14. Done in the refactor + known limitations

Done:
- `pipeline/` (config/stages/runner/readiness/cli) + `run_automag.py` +
  `input_automag.py`; single command, arbitrary directory, results in the
  launch dir.
- **Stage instances**: a stage may be repeated with different calculators;
  per-instance work dirs, `coll_calculator`, `depends_on`, and the
  `reuse`/`on_missing`/`on_unsupported` policies.
- Calculator-aware `2_coll` lookups (`common/results.coll_bundle`) so a stage
  can read another calculator's results.
- `2_relax` stage (most stable `2_coll` configuration, converged moments).
- `calculators/` registry + VASP/QE/FHI-aims/prophet/generic interfaces;
  capability gate skips unsupported instances.
- `--wait` DAG scheduler over `(instance, phase)` with file-based readiness.
- Automatic VAMPIRE/ESpinS run (`3_monte_carlo/5_run_mc.py`).
- Path decoupling (`common/paths.py`), generic result writer, copyable
  examples incl. `examples/mixed_*`.
- Ninja/PyInstaller single-file build with `--run-script`/`--run-helper`.
- `input.{py,yaml,json}` for the controlling file and the stage settings.
- Optional `spin_values`/`magnetic_atoms` (pymatgen `DEFAULT_MAGMOMS`) and a
  second enumeration backend (`enumeration_mode = "pymatgen"`, i.e.
  `MagneticStructureEnumerator`).
- External codes (enumlib/VAMPIRE/ESpinS) built by `./build.sh --external` and
  registered in `external/codes.json`.
- **Execution layer**: hosts (`local`/`ssh`/`shell`, per-instance `host`),
  smart launcher selection (`sbatch`/`qsub`/process from the `header` or the
  host scheduler), `--once`/`--wait`/`--daemon`/`--status`/`--report`, SQLite
  state with resuming, parallel calculation jobs (`max_parallel`), retries
  (`max_attempts`), and re-running only the steps whose settings changed
  (`config_hash` + `reset_step`).
- **Parameter sweeps** in `mae` (`params` values given as lists) with the
  `MAE_convergence_*` / `MAE_curve_convergence_*` summaries (`common/sweep.py`).

Known limitations / possible next work:
- `lin_response`, `relax`, `mae` are VASP-only (`relax_params` are VASP tags,
  MAE needs non-collinear support). Other backends need parameter mapping or
  new stage implementations.
- `SubmitFirework` path is stale.
- QE AFM per-species limitation (above).
- Suffix inconsistency in `write_charges.py`.
- A step whose calculations are missing results is marked `done` with a warning
  rather than failed (deliberate: a non-converged magnetic configuration must
  not stop the pipeline); the dependent analysis step then fails with the real
  diagnostic.
- The sweep expansion covers `4_mae` only; `conv_tests` has its own value lists
  (`encut_values`/`sigma_values`/`kpts_values`).
- Two instances of the same stage with the same calculator still share
  `CalcFold` files (results are appended twice).
