## Python environment

All commands in terminal must be run as:

    source /home/dpoletaev/miniconda3/bin/activate automag && <command>

or:

    source /home/dpoletaev/miniconda3/etc/profile.d/conda.sh && conda activate automag && <command>

Never use the system Python or any other environment.
Conda binary: /home/dpoletaev/miniconda3/bin/conda.

## automag agent quickstart

Full reference: `docs/AGENT_GUIDE.md` — open it only for pipeline/calculator/stage work.
Session state, decisions and known pitfalls: `docs/SESSION_HANDOFF.md` — read §0 and
§10–§12 before changing anything.

- Entry: `run_automag.py`; config: `input_automag.{py,yaml,json}`
  (`PIPELINE` dict; priority `.py` → `.yaml` → `.json`). Stage settings load
  the same way (`input.{py,yaml,json}` via `common.input_loader`).
- No config file in the launch dir → `pipeline/bootstrap.py` writes a default
  `input_automag.py` when a structure file is present (plan: `coll` with the
  bundled prophet_spin; `4_mae` is deliberately excluded — the released model is
  isotropic and would give MAE = 0).
- Run all: `cd <project dir> && python <automag>/run_automag.py --wait`;
  plan only: add `--dry-run`. Results (stage dirs, `CalcFold`) land in the
  launch dir; never write inside the automag checkout.
- A stage may repeat with different calculators (instances in `PIPELINE["run"]`:
  `id`, `stage`, `calculator`, `coll_calculator`, `reuse`, `on_missing`).
- Never rename stage dirs or `CalcFold` result filenames; post-processing depends on them.
- Read `2_coll` results only via `common.results` (`coll_bundle`,
  `find_state`, `most_stable_state`) and always pass the calculator.
- `spin_values`/`magnetic_atoms` are optional (pymatgen `DEFAULT_MAGMOMS` via
  `common.magnetic`); `2_coll` also supports `enumeration_mode="pymatgen"`.
- New calculator: add `calculators/<name>.py` + register in `calculators/__init__.py`;
  do not edit pipeline/stage code.
- Use `common.paths.code_root()` / `data_root()`; `AUTOMAG_PATH` is legacy.
  Prophet checkpoints: `common.paths.models_dir()` / `model_file(name)` —
  `models/prophet-v1-mag.pt` and `models/prophet-v1-mag-head.pt` are in the
  repo and also inside the frozen bundle.  Never hard-code a model path in an
  example: leave `model_path` out (the interface falls back to the bundle).
- Validate: `python -m py_compile ...`, `python tools/smoke_test.py`,
  `python tools/execution_test.py`, `python tools/mae_chain_test.py`,
  `python tools/launchers_test.sh` (all without DFT), and
  `cd examples/vasp && python ../../run_automag.py --dry-run --wait`
  (plus `examples/prophet`, `examples/prophet_then_vasp`).
- Build: `./build.sh [--onedir|--all|--external]` (ninja + PyInstaller) →
  `dist/run_automag`.  `prophet_spin` + torch + `models/` are bundled by
  default (~5.7 GB); `--no-prophet` gives the ~490 MB bundle.  The frozen build
  must keep working without conda (`env -i PATH=/usr/bin:/bin <bin> --dry-run`).
- External codes (enumlib/VAMPIRE/ESpinS): `./build.sh --external
  [--<code>-path P]` writes `external/codes.json`; use `common.external`, never
  hard-coded paths.
