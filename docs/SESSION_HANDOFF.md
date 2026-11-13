# automag — session handoff

State of the repository after the multi-calculator / instances / inputs / external-codes
/ enumeration session. Read §0, then only the section you need.
Code is the source of truth; this file records *why* things are the way they are.

---

## 0. TL;DR for the next agent

- Everything listed in §13 "verified" was executed successfully at the end of the
  session. Nothing was left half-applied (the tree compiles).
- **Before changing anything** run §12 (3 commands, ~1 min) to confirm the baseline.
- The most dangerous areas are §10 (bugs already found — do not regress) and §11
  (tooling pitfalls that wasted a lot of time).
- Two items were **never executed** because the tools are absent in this env:
  `./build.sh` (needs `ninja` + `PyInstaller`) and real DFT runs. See §13/§14.
- Git state is **uncommitted** (§15). Do not assume `git pull`-clean.

---

## 1. Environment facts

| Fact | Value |
| --- | --- |
| Python env (mandatory for every command) | `source /home/dpoletaev/miniconda3/bin/activate automag && <cmd>` |
| Conda binary | `/home/dpoletaev/miniconda3/bin/conda` |
| Repo root | `/home/dpoletaev/magnetic/automag` |
| Toolchain present | `git`, `make`, `cmake`, `gfortran`, `mpirun` |
| **Absent** | `ninja`, `PyInstaller`, Intel Fortran/MKL |
| Python deps present | pymatgen, ase, `PyYAML 6.0.2`, numpy, matplotlib, spglib |
| Network | works (GitHub clones succeeded) |
| `terminal` tool noise | prints `CondaError: Run 'conda init'` + `screen size is bogus` — **harmless**, commands run |

`terminal` tool constraints that cost time:
`cd` must point inside the project root (no `/tmp/.../repo`); no `$VAR` inside the
command (the shell mangles it); `grep -c` exits 1 on zero matches and breaks `&&`.

---

## 2. Repo map — what is new or heavily changed in this session

| Path | Status | Role |
| --- | --- | --- |
| `pipeline/runner.py` | **rewritten** | instance scheduling, policies, auto-PYTHONPATH, frozen `--run-script` |
| `pipeline/config.py` | rewritten | `StageInstance`, work dirs, `external_codes`, `_check_workdir`, config path resolution |
| `pipeline/readiness.py` | rewritten | `stage_work_dir`, `coll_calculator`, `"lowest"` for MAE |
| `pipeline/stages.py` | changed | `work_dir` names; `spin_values` dropped from `coll.required` |
| `pipeline/cli.py` | changed | `--instances`, config default `None`, warns about unavailable codes |
| `common/results.py` | rewritten | `coll_bundle`, `available_coll_calculators`, calculator-aware lookups, masked spin patterns |
| `common/input_loader.py` | rewritten | `.py` → `.yaml`/`.yml` → `.json`, `apply_stage_input` |
| `common/external.py` | **new** | registry of enumlib/VAMPIRE/ESpinS (+ pre/post commands) |
| `common/magnetic.py` | **new** | DEFAULT_MAGMOMS resolution + `MagneticStructureEnumerator` backend |
| `common/paths.py` | rewritten | `frozen()`, `code_root()` = `_MEIPASS` when frozen, `helper_command()` |
| `common/SubmitManual.py` | changed | auto-submit of `run_all_*.sh`, frozen helper calls |
| `run_automag.py` | rewritten | `--run-script` / `--run-helper` |
| `2_relax/` | **new** | re-relaxation stage (`1_submit.py`, `2_collect.py`, `input_template.py`) |
| `3_monte_carlo/5_run_mc.py` | **new** | runs VAMPIRE/ESpinS from the config |
| `0_conv_tests`, `1_lin_response`, `2_coll`, `2_relax`, `3_monte_carlo`, `4_mae` scripts | changed | `apply_stage_input`, `coll_calculator`, `common.external`, `common.magnetic` |
| `build.ninja`, `build.sh`, `tools/` | **new** | ninja + PyInstaller build, external-code build, smoke test |
| `examples/` | **new/changed** | 9 copyable examples, unified `run_pipeline.sh` |
| `docs/AGENT_GUIDE.md` | updated | the dense reference (§0–§14) |
| `docs/REPO_MAP.md` | **deleted** | overlapped with AGENT_GUIDE; do not re-create |

Pre-existing (NOT from this session, do not "clean up"): legacy per-project dirs inside
stage folders (`2_coll/Al10Fe20Mo10`, `3_monte_carlo/Fe12O18`, `4_mae/Fe4O6_alpha`, …)
and the already-modified `3_monte_carlo/Fe12O18/input.py`, `3_monte_carlo/example_input.py`,
`common/write_charges.py`, `common/write_output.py`, `common/write_output_vasp.py`.

---

## 3. Execution model (instances)

```
input_automag.{py,yaml,json}  ->  PipelineConfig.instances[StageInstance]
                                     |
                                PipelineRunner
    non-wait: each instance -> its phases (scripts, in order)
    --wait:   DAG over (instance_id, phase) -> run -> poll(readiness) -> continue
                                     |
        <workdir = launch dir by default>/
            0_conv/ 1_lin_response/ 2_coll/ 2_relax/ 3_monte_carlo/ 4_mae/
            CalcFold/  geometries/
```

- `Stage.work_dir` names the sub-directory; the **code** directories keep their old
  names (`2_coll` code ↔ `2_coll` work dir; `conv_tests` code ↔ `0_conv` work dir).
- One instance of a stage → `<workdir>/<Stage.work_dir>`; several → `…/<instance id>`.
- `CalcFold` and `geometries` are always global to the workdir.
- `workdir` defaults to `Path.cwd()`; `config._check_workdir` rejects the automag root
  and anything inside the checkout that is not under `examples/`.

Runner internals worth knowing:
`Task(instance, phase, scripts, depends, wait, wait_before)`, key `(instance.id, phase)`.
Cross-instance dependency = dependency instance's **last planned phase**.
`CAPABILITIES = {lin_response: supports_ldau, relax: supports_relaxation, mae: supports_noncollinear}`.
`_stage_done()` markers: coll → `most_stable_state`; relax → `relaxed_*.json`;
monte_carlo → non-empty `output`; mae → `MAE_curve_*.txt`; conv → `CalcFold/<f>_<mode>_<calc>.txt`;
lin → charges file.

Policies: `reuse` (default on for dict `monte_carlo`/`mae`), `on_missing`
(`skip`/`error`), `on_unsupported` (`skip` default for dict instances, `error` for legacy).

---

## 4. Config formats and merge order

Controlling file: `input_automag.py` → `.yaml`/`.yml` → `.json` (first that exists).
Stage settings: `input.py` → `input.yaml` → `input.json` (`common.input_loader.resolve_input`).
YAML/JSON may wrap in `PIPELINE:`/`pipeline:` or be flat.

Merge order for the generated `input.py`:
`submit`/`common` → `calculators[<calc>]` → `stages[<stage>]` → instance dict,
then always `calculator`, `poscar_file`, `struct_suffix`, `coll_calculator`.

`submit`/`common` keys: `calculator`, `use_fireworks`, `jobheader`, `calculator_command`,
`environment_activate`, `environment_deactivate`, `parallel_over_configurations`,
`auto_submit`, `exec_command`, `external_codes`.

Instance keys: `id`, `stage`, `calculator`, `work_dir`, `depends_on`, `coll_calculator`,
`coll_instance`, `on_missing`, `on_unsupported`, `reuse`, `enabled`, `steps`, `description`;
everything else is a stage setting.

YAML/JSON cannot express `range(...)` — use a list; stage code only iterates.

---

## 5. Data contracts (do not change)

`CalcFold/<formula><suffix>_singlepoint_<calc>.txt` — one header line per configuration,
next line `initial_magmoms=[…]  final_magmoms=[…]`. `values[-4]` is parsed by
`2_coll/2_plot_results.py`; `final_magmoms=` is parsed by `common.results.final_magmoms`.

`2_coll` bundle (search recursively — one per instance):

```
2_coll[/<instance>]/<formula><suffix>_<calc>/
    trials_<formula><suffix>/{settingNNN.vasp, configurationsNNN.txt}
    <formula><suffix>_<calc>_{states,energies,setting}NNN.txt/.vasp
    <formula><suffix>_<calc>_most_stable_configuration.txt
    <formula><suffix>_<calc>_configurations_stability_ascending.txt
```

`configurationsNNN.txt` lines: `<state>  <spin per magnetic atom>`, states
`fm<N>`/`afm<N>`/`fim<N>`/`nm` — `4_mae` and MC depend on that naming.

Other: `3_monte_carlo/output` (2 columns T,M), `MAE_theta_phi_*.txt`,
`MAE_curve_*.txt`, `settingNNN_<conf>_processed.{vasp,json}`.

Always read `2_coll` through `common.results` **with an explicit calculator**:
`coll_bundle(root, formula, suffix, calculator)`, `most_stable_state(...)`,
`find_trials_dir(...)`, `find_state(..., calculator=…)`, `load_state(..., calculator=…)`.

---

## 6. Cross-calculator rules (the core feature)

- `coll_calculator` = whose `2_coll` bundle an instance reads; default = its own
  `calculator`. Downstream stages must use it for the bundle, for `find_state`, for
  `CalcFold/..._singlepoint_<coll_calculator>.txt`, and for `"lowest"` resolution.
- `"lowest"` in `monte_carlo`/`mae`/`relax` resolves per `coll_calculator`
  (`common.results.resolve_configuration`).
- Unsupported calculator for a stage → skipped with a warning; the runner additionally
  prints `note: 2_coll results exist for calculator(s): …` and suggests declaring the
  missing instance (`runner._warn_alternatives`).
- An instance whose `coll_calculator` has no declared coll instance falls back to
  checking data on disk (`runner._coll_data_available`) and is skipped/errors per
  `on_missing`.
- Mixed-and-repeated pipelines are the primary use case; see
  `examples/mixed_prophet_vasp` (prophet → VASP → prophet → VASP).

---

## 7. Magnetic moments and enumeration (latest feature)

`common/magnetic.py`:

- `resolve_spin_values(structure, spin_values, magnetic_atoms, overrides=None)`
  → `{symbol: [magnitude, ...]}`; missing entries filled from
  `pymatgen.analysis.magnetism.analyzer.DEFAULT_MAGMOMS`; untabulated elements are
  omitted ⇒ non-magnetic. **Integral moments stay `int`** (stage scripts format them
  with `{:2d}` and use them as spin multiplicities).
- `resolve_magnetic_elements(...)` — explicit `magnetic_atoms` wins, else non-zero moment.
- `structure_moments(...)`, `element_symbol(site)`, `site_spin(site)`.
- `ensure_enumlib_on_path()` — pymatgen's enumerator looks for `enum.x`/`makestr.x`
  on `PATH`; this adds the directories from `external/codes.json`.
- `pymatgen_configurations(...)` → `(lattices, coordinates, configurations, origins)`.

`2_coll/1_submit.py` selects the backend:

| `enumeration_mode` | behaviour |
| --- | --- |
| `"enumlib"` (default) | original Wyckoff driver + optional sublattice split + `enum.x`/`makeStr.py` |
| `"pymatgen"` | `MagneticStructureEnumerator` with `enumeration_strategies` and `max_orderings` |

Both fill the same three lists and write the same `trials_*/` layout, so nothing
downstream changed. `enumlib` writes multiplicities (`5 -5 …`), `pymatgen` writes signs
(`1 -1 …`) — both are consumed via `np.sign`, keep it that way.

Strategies: `ferromagnetic`, `antiferromagnetic`, `antiferromagnetic_by_motif`,
`ferrimagnetic_by_motif`, `ferrimagnetic_by_species`, `nonmagnetic`.

---

## 8. External codes (enumlib / VAMPIRE / ESpinS)

Runtime registry `common/external.py`; resolution order (later wins):
built-in defaults (`enum.x`/`makeStr.py`, `vampire`, `mc.x` on `PATH`) →
`code_root()/external/codes.json` → `<exe dir>/external/codes.json` (frozen) →
`AUTOMAG_EXTERNAL_CODES` (file or dir) → `external_codes` in the pipeline input.

Entry: `{available, command | enum+makestr+python, pre_command, post_command}`.
`pre_command`/`post_command` run in one bash around the executable (so `module load`,
`source`, `conda activate` work). `.py` executables are launched with the recorded
`python` interpreter, not their shebang.

Build/register: `python tools/build_external.py` = `./build.sh --external`, ninja target
`external`. `--<code>-path P` (+ `--<code>-pre-command`/`--<code>-post-command`) or
`--config FILE` (JSON/YAML) records an existing installation and **skips** cloning and
compiling. Each build step is retried once with `system_env()` (conda `g++` fails to link
against the system glibc). Failure ⇒ `available: false` + warning telling the user to
compile manually and pass the path; the rest of the pipeline keeps working, only the
dependent step stops with the same instruction.

Current state of `external/codes.json` on this machine (git-ignored):

| code | state |
| --- | --- |
| enumlib | built, `external/src/enumlib/src/enum.x` |
| vampire | built, `external/src/vampire/vampire-serial` |
| espins | `available: false` — ESpinS needs Intel Fortran + MKL (`make.sys`) |

Recipes: enumlib = `git clone --recursive` + `make F90=gfortran` in `symlib/src`, then
`make enum.x makestr.x F90=gfortran`; VAMPIRE = `make serial` → `vampire-serial`;
ESpinS = `make` → `mc.x`.

---

## 9. Build (ninja, PyInstaller, frozen)

- `./build.sh` → ninja → `dist/run_automag`; `./build.sh --external [opts]`;
  `./build.sh --all`; `./build.sh --onedir`; `./build.sh --clean`.
- `build.ninja` includes the generated `build/deps.ninja`
  (`tools/write_ninja_deps.py --output …`) because ninja has no globbing;
  `builder` target is `generator = 1`.
- `tools/build_executable.py` = PyInstaller `--onefile` with `--add-data` for all stage
  dirs + `common`/`calculators`/`pipeline`/`geometries`.
- Frozen runtime: `code_root()` = `sys._MEIPASS`; the runner starts stages as
  `<exe> --run-script <dir>/<script>`; job scripts and the prophet interface use
  `common.paths.helper_command()` → `<exe> --run-helper common/<script>`.
  Paths are **relative to the bundle root** because every launch of a one-file build
  unpacks to a fresh temp dir.

---

## 10. Bugs found in this session — do not regress

1. **`from input import *` needs cwd on `sys.path`.** `python /abs/script.py` puts the
   *script's* dir on `sys.path`, so the runner adds `inst.work_dir` to `PYTHONPATH`
   (`runner._run_script`). Symptom if broken: `ModuleNotFoundError: No module named 'input'`.
2. **Nothing executed the generated `run_all_*.sh`.** `SubmitManual.submit()` now flushes
   them (`sbatch`/`qsub`/`sh` from the job header), idempotently: one process-level
   buffer per script, only commands queued since the last flush are launched, the file is
   truncated on first touch. Disable with `"submit": {"auto_submit": false}`
   (runner exports `AUTOMAG_AUTO_SUBMIT=0`). Requires `parallel_over_configurations=True`,
   otherwise nothing is submitted (a warning is printed).
3. **MC/MAE instance without a declared coll instance ran anyway** and crashed inside the
   script. Fixed via `_coll_data_available` + `on_missing`.
4. **`site.spin` does not exist** on `PeriodicSite`; the spin is `site.specie.spin`
   (encoded as `Fe,spin=np.float64(-5.0)`).
5. **`site.specie` may be `Element`, which has no `.element`** → use
   `magnetic.element_symbol(site)`.
6. **Float moments broke `{:2d}`** — keep integral magnitudes integral (`_numeric`).
7. **pymatgen groups configurations by lattice, not by lattice+frac_coords**, otherwise
   AFM patterns of one cell land in a second `setting`.
8. **`--wait` + `auto_submit=false` hangs forever** (it polls for calculations that never
   start). Use non-wait when testing, or a small `wait_timeout`.

---

## 11. Tooling pitfalls (learned the hard way)

- Large `write_file` calls were occasionally mangled/not applied. Prefer several small
  files/edits and verify each with `grep -n` or `python -m py_compile`.
- Editing a block that starts with `if …:` — verify the indentation survived
  (`sed -n 'A,Bp'`) before moving on; one edit once injected an import into the wrong
  branch and another replaced a `def` line.
- For multi-file mechanical edits, write a throw-away script, run it, `py_compile`,
  `grep` the result, then delete the script. That is how the 13 stage prologues were
  migrated to `common.input_loader.apply_stage_input`.
- `git --no-pager` for read-only git; `GIT_EDITOR=true` for anything that may open an editor.
- `git check-ignore -v <path>` to confirm generated dirs stay out of the repo.

---

## 12. Validation protocol (run before and after any change)

```bash
cd /home/dpoletaev/magnetic/automag
source /home/dpoletaev/miniconda3/bin/activate automag && python -m py_compile \
  run_automag.py input_automag.py pipeline/*.py calculators/*.py common/*.py \
  0_conv_tests/*.py 1_lin_response/*.py 2_coll/*.py 2_relax/*.py \
  3_monte_carlo/1_coupling_constants.py 3_monte_carlo/2_write_vampire_ucf.py \
  3_monte_carlo/4_write_espins_mcin.py 3_monte_carlo/5_run_mc.py \
  4_mae/1_submit.py 4_mae/2_plot_results.py 4_mae/3_submit.py 4_mae/4_plot_results.py \
  tools/*.py examples/*/input_automag.py && echo COMPILE_OK

source /home/dpoletaev/miniconda3/bin/activate automag && python tools/smoke_test.py
#   -> "smoke test OK: cross-calculator reads, reuse/skip policies, external-code wiring"

source /home/dpoletaev/miniconda3/bin/activate automag && \
# (examples are not frozen; list them from the directory when they change)
for d in vasp qe fhi-aims prophet mixed_fhi_aims_vasp mixed_qe_vasp \
         pymatgen_enumeration remote_vasp remote_fhi_aims remote_qe; do
  (cd examples/$d && python ../../run_automag.py --dry-run --wait >/dev/null && echo "OK $d")
done
```

`tools/smoke_test.py` needs no DFT and no enumlib: it fabricates a `prophet_spin`
`2_coll` bundle in a temp dir and asserts real script execution, cross-calculator reads,
`reuse`, `on_missing`, `auto_submit=False`, and external-code wiring (fake VAMPIRE with
pre/post commands, plus the `available: false` error path).

Clean up after testing — these are git-ignored but must not linger in `examples/`:
`examples/*/{0_conv,1_lin_response,2_coll,2_relax,3_monte_carlo,4_mae,CalcFold,geometries,__pycache__}`, `runs/`.

---

## 13. Verified at end of session

- `py_compile` of the whole tree (incl. examples' configs) — OK.
- `tools/smoke_test.py` — OK (see above).
- `--dry-run --wait` on all 9 examples — OK, including the mixed sequences.
- Mixed plan for the user's target sequence
  `2_coll(prophet) → relax(vasp) → 2_coll(vasp) → mae(prophet) → mae(vasp) → mc(prophet) → mc(vasp)`:
  `mae_prophet` skipped (no non-collinear support) with a warning; `mc_vasp` skipped when
  VASP data is absent, with the advisory note.
- enumlib: cloned + compiled with gfortran, `2_coll/1_submit.py` really invoked
  `enum.x` + `makeStr.py` and produced `trials_*/enumlib1/`.
- VAMPIRE: compiled after the automatic system-toolchain retry; `5_run_mc.py` launched a
  (fake and real-path) VAMPIRE through the registry, pre/post commands executed,
  `output` copied.
- ESpinS: build fails as expected (no MKL) → `available: false` + warning + actionable
  runtime error mentioning `./build.sh --external`.
- `input_automag.{py,yaml,json}` and stage `input.{yaml,json}` all load; priority
  `.py` > `.yaml` > `.json` verified.
- `common.magnetic` defaults: Fe4O6 → `{'Fe': [5]}`, `['Fe']`; explicit
  `{'Fe': [4]}` → `{'Fe': [4]}`; non-tabulated element stays non-magnetic.
- pymatgen backend: Fe4O6 → `fm1`, `afm1..afm3` in one `setting001.vasp`;
  `enumeration_strategies` + `max_orderings` honoured.

## Not verified (tools absent / out of scope)

- `./build.sh` and `./build.sh --all` end-to-end (no `ninja`, no `PyInstaller`).
  Only `tools/build_executable.py --dry-run` and `tools/write_ninja_deps.py` were run.
  Therefore `dist/run_automag` and the frozen `--run-script`/`--run-helper` paths are
  **code-reviewed but not executed**.
- Any real VASP / QE / FHI-aims / prophet calculation.
- `2_coll` end-to-end with the VASP calculator (no `POTCAR`/`VASP_PP_PATH` here) — the
  collinear generation itself was validated with `calculator="prophet_spin"`.
- `3_monte_carlo/4_write_espins_mcin.py` execution (needs ESpinS inputs) — only its
  command-resolution path was reviewed.

---

## 14. Known limitations / sensible next steps

1. `lin_response`, `relax`, `mae` are VASP-only in practice (`relax_params` are VASP tags,
   MAE needs non-collinear support). Doing them for QE/FHI-aims means mapping parameters
   and is a real feature, not a fix.
2. Two instances of the same stage **with the same calculator** still share `CalcFold`
   result files (results appended twice). A per-instance suffix in the result file name
   would fix it but changes §5 contracts.
3. `SubmitFirework` is stale (stage scripts pass `struct_suffix` it does not accept);
   `use_fireworks=False` is the supported path.
4. `write_charges.py` writes `<formula>_charges_<calc>.txt` without `struct_suffix` while
   `1_lin_response/2_plot_results.py` reads it with the suffix — consistent only for an
   empty suffix.
5. QE `starting_magnetization` is per species; AFM within one species needs distinct
   sublattice labels.
6. `--wait` has no scheduler integration (polls result files); file-count checks continue
   after `wait_stable_seconds` without progress.
7. ESpinS stays unavailable without Intel Fortran + MKL — document the manual path.

---

## 15. Git state (uncommitted)

Tracked files modified: `.gitignore`, `0_conv_tests/{1_submit,2_plot_results}.py`,
`1_lin_response/{1_submit,2_plot_results}.py`, `2_coll/{1_submit,2_plot_results}.py`,
`3_monte_carlo/{1_coupling_constants,2_write_vampire_ucf,4_write_espins_mcin}.py`,
`4_mae/{1_submit,2_plot_results,3_submit,4_plot_results,input_template}.py`,
`2_coll/input_template.py`, `2_relax/input_template.py`, `AGENTS.md`, `README.md`,
`common/SubmitManual.py`.

Untracked (new): `2_relax/`, `3_monte_carlo/5_run_mc.py`, `build.ninja`, `build.sh`,
`calculators/`*, `common/{external,input_loader,magnetic,paths,results,run_prophet,write_output_ase}.py`,
`docs/`, `examples/`, `input_automag.py`, `pipeline/`, `run_automag.py`, `tools/`.

(\* `calculators/` and some `common/` files were created in earlier turns of the same
work and are still untracked.)

Ignored: `runs/`, `/build/`, `/dist/`, `/external/`, `*.spec`, `examples/*/<stage dirs>`.

Nothing is committed; the next agent should commit or review with `git diff` first.
