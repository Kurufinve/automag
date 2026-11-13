# Automag

An automatic workflow for magnetic materials: it finds the ground collinear
magnetic state of a structure, computes the Heisenberg exchange couplings and
estimates the critical temperature of the magnetically ordered → paramagnetic
transition, and evaluates the magnetocrystalline anisotropy energy (MAE).

Everything is driven from **one input file and one command**.  Automag decides
what has to run next, starts the calculations itself (in-process, as separate
processes, or through a batch scheduler on this or another machine), stores the
progress in SQLite, and can be stopped and resumed at any point.

```
structure ──► conv_tests ──► coll ──► relax ──► monte_carlo ──► Tc
                   │           │                      ▲
                   └─ lin_response (U)                │
                               └──────────────────► mae ──► MAE, K1, K2
```

| stage | directory | what it does |
| --- | --- | --- |
| `conv_tests` | `0_conv_tests` | ENCUT / k-mesh / SIGMA convergence tests |
| `lin_response` | `1_lin_response` | Hubbard U by linear response |
| `coll` | `2_coll` | enumeration and single-point energies of collinear magnetic configurations |
| `relax` | `2_relax` | re-relaxation of the most stable configuration |
| `monte_carlo` | `3_monte_carlo` | exchange couplings (Heisenberg model) + VAMPIRE/ESpinS run |
| `mae` | `4_mae` | magnetocrystalline anisotropy: θ/φ scan, easy/hard axis, K1, K2 |

## What this version can do

* **One controlling file** (`input_automag.py|yaml|json`) for every stage, with
  per-instance overrides and several calculators in one run.
* **Three ways to run it**: supervised (`--wait`/`--daemon`), one-shot
  (`--once`, cron-friendly) and step-by-step (`--steps`, the historical
  behaviour).
* **Resumable**: the state of every step and every started job lives in
  `automag_state.sqlite`; an interrupted run continues where it stopped, and a
  changed setting re-runs only the steps it affects.
* **Anywhere**: a step can run in the automag process, as a separate local
  process, as a local batch job (SLURM/PBS), or on a remote machine (also via
  its scheduler).  The choice is per step in the input file.
* **Parallel**: the independent calculations of a step are started as separate
  jobs with a configurable limit; `4_mae` respects the `z → θ/φ → curve` order
  on its own.
* **Convergence studies**: give a parameter of `4_mae` as a list
  (`"encut": [400, 500, 600]`) and every combination becomes an independent MAE
  calculation, compared in a convergence summary.
* **Calculator-agnostic**: VASP (reference), Quantum ESPRESSO, FHI-aims, ASE
  calculators and custom plugins; a stage may read the `2_coll` results of a
  *different* calculator.
* **Works out of the box**: the Prophet-V1-Mag checkpoints ship in `models/`,
  so a directory with nothing but a structure file gets a working DFT-free plan
  (`2_coll` with the neural-network model) without any input file.
* **Self-contained binary**: the compiled build carries Python, the stage
  scripts, the model checkpoints and `prophet_spin`, and runs on a machine
  without conda.

Installation and the optional single-file build: **[docs/INSTALL.md](docs/INSTALL.md)**.
Developer/agent reference: **[docs/AGENT_GUIDE.md](docs/AGENT_GUIDE.md)**.

## Installation

Python **3.10+** is required (3.12 is tested).  Pick one of the three ways; all
three are described in detail in [docs/INSTALL.md](docs/INSTALL.md).

### 1. conda (recommended)

```bash
conda create -n automag python=3.12 -y
conda activate automag
git clone https://github.com/Kurufinve/automag.git && cd automag

# DFT calculators, stage scripts, pipeline
python -m pip install -r requirements.txt

# optional: the Prophet calculator (torch, ~2 GB download)
python -m pip install -r requirements-prophet.txt
```

### 2. venv without conda

```bash
python3 -m venv ~/venvs/automag && source ~/venvs/automag/bin/activate
cd automag
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install -r requirements-prophet.txt      # optional
```

### 3. System Python (no conda, no venv)

```bash
# RHEL/Rocky: sudo dnf install python3.12 python3.12-devel gcc gcc-c++ #                          blas-devel lapack-devel
# Debian/Ubuntu: sudo apt install python3 python3-venv python3-dev g++ #                          libblas-dev liblapack-dev
cd automag
python3 -m pip install --user -r requirements.txt
python3 -m pip install --user -r requirements-prophet.txt   # optional
```

`enumlib` (the magnetic-configuration enumerator used by `2_coll`) and the
Monte Carlo backends are external programs, not Python packages:

```bash
./build.sh --external            # download + compile enumlib / VAMPIRE / ESpinS
./build.sh --external --espins-path /path/to/mc.x   # or register existing ones
```

Without them `2_coll` stops with a clear message; see
[INSTALL §6](docs/INSTALL.md).

### Check the installation

```bash
python -m py_compile run_automag.py                # all sources parse
python tools/smoke_test.py                        # fast, no DFT code needed
cd examples/vasp && python ../../run_automag.py --dry-run --wait
```

## Run from source

```bash
source /path/to/miniconda3/bin/activate automag     # or: source ~/venvs/automag/bin/activate
cd my_project                                       # input_automag.py + structure.vasp
python /path/to/automag/run_automag.py --dry-run --wait   # show the plan, run nothing
python /path/to/automag/run_automag.py --wait             # run everything and wait
```

Results are written next to `input_automag.py`: the stage directories
(`0_conv/`, `2_coll/`, `2_relax/`, `3_monte_carlo/`, `4_mae/`),
`CalcFold/`, `geometries/` and `automag_state.sqlite`.  The automag checkout is
never written to.

### No input file at all

```bash
mkdir Fe4O6 && cd Fe4O6
cp /somewhere/structure.vasp .
python /path/to/automag/run_automag.py --wait
```

If the directory holds nothing but a structure (`structure.vasp`, `POSCAR`,
`*.vasp`, `*.cif`, `*.poscar`), automag writes a default `input_automag.py` and
continues with it: the magnetic search (`2_coll`) evaluated with the bundled
Prophet-V1-Mag model - no DFT code, no cluster, no input file.  The generated
file is an ordinary input file; edit it to change anything (a `2_coll` of this
size takes a few minutes per configuration on CPU).  `4_mae` is **not** part of
the default plan: the released Prophet weights are isotropic, so the
anisotropy would be zero (see `examples/prophet_mae/` if you need the stage).

### Run modes

| command | behaviour |
| --- | --- |
| `--wait` | run the whole pipeline, waiting for the calculations between steps |
| `--daemon` | same, but keep running and supervise |
| `--once` | advance until a step must wait, record it, exit (cron-friendly) |
| `--status` / `--report` | print the stored state / a summary and exit |
| `--dry-run` | print the plan, execute nothing |
| `--steps submit,postprocess`, `--stages coll,mae`, `--instances coll_vasp` | restrict what runs |
| `--execution process`, `--fresh`, `--state-file P`, `--poll-interval S` | override execution/state settings |
| `--list-stages`, `--list-hosts` | what the pipeline knows about |

`--wait`/`--once` resume from `automag_state.sqlite`: an interrupted run
continues where it stopped, and a changed setting re-runs only the steps it
affects.

### Launchers shipped with every example

`examples/*/` contain ready-made scripts: copy the directory anywhere, edit
`input_automag.py`, and run one of

```bash
./run_pipeline.sh                  # sources, one command, waits until finished
./run_pipeline_once.sh             # sources, one pass (manual restarts, cron)
./run_pipeline_once_loop.sh        # the same passes in a loop
./run_pipeline_compiled.sh         # compiled binary (AUTOMAG_EXE)
./run_pipeline_slurm.sh            # sources, the supervisor itself as a SLURM job
./run_pipeline_slurm_compiled.sh   # compiled binary, submitted to SLURM
```

See [examples/README.md](examples/README.md) for the environment variables
(`AUTOMAG_HOME`, `AUTOMAG_EXE`, `AUTOMAG_CONDA_ENV`, ...).

## The controlling file

`input_automag.py` in the launch directory (or `input_automag.yaml` /
`input_automag.json` — the order is `.py` → `.yaml`/`.yml` → `.json`; the same
fallback applies to the per-stage `input.*` files).  Minimal example:

```python
VASP_PARAMS = {
    "xc": "PBE", "prec": "Accurate", "ncore": 4, "encut": 600,
    "ediff": 1e-6, "ismear": 1, "sigma": 0.1, "nelm": 200, "kpts": 20,
    "lcharg": False, "lwave": False, "isym": 0,
    "ldau": True, "ldautype": 2, "ldaul": [2, -1],
    "ldauu": [4.0, 0], "ldauj": [0, 0], "ldauprint": 2,
}

PIPELINE = {
    "structure": "structure.vasp",
    "run": ["conv_tests", "coll", "relax", "monte_carlo", "mae"],

    "submit": {
        "calculator": "vasp",
        "use_fireworks": False,
        "jobheader": """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=YOUR_PARTITION
""",
        "calculator_command": "mpirun vasp_std",
        "environment_activate": "source /path/to/conda.sh && conda activate automag",
    },

    "coll": {"supercell_size": 1, "spin_values": {"Fe": [4]}, "params": VASP_PARAMS},
    "relax": {"params": VASP_PARAMS, "relax_params": {"ibrion": 2, "isif": 3, "nsw": 100}},
    "monte_carlo": {"configuration": "lowest", "cutoff_radius": 4.1},
    "mae": {"configuration": "lowest", "Nph": 8, "Nth": 4, "N_MAE": 8,
            "params": VASP_PARAMS},
}
```

A fully commented version, including `conv_tests`, is
`examples/vasp/input_automag.py`; `examples/README.md` lists the ready-made
setups.

### Top-level keys

| key | meaning |
| --- | --- |
| `structure` | **required**; absolute, relative to the input file, or a name in `automag/geometries` |
| `workdir` | where results are written (default: the launch directory) |
| `struct_suffix` | suffix added to the `CalcFold` compound directory |
| `run` | stages to execute, in order — names or instance dictionaries |
| `submit`, `common` | settings shared by every stage (calculator, job header, environment, `use_fireworks`, `auto_submit`, `parallel_over_configurations`, `external_codes`) |
| `hosts` | named execution hosts (see [Where a step runs](#where-a-step-runs)) |
| `execution` | default execution plan: mode, `max_parallel`, `max_attempts`, per-phase overrides |
| `max_parallel_tasks`, `max_attempts`, `poll_interval` | short forms of the above |
| `state_file` | SQLite file (default `<workdir>/automag_state.sqlite`) |
| `wait_timeout`, `wait_interval`, `wait_stable_seconds` | bounds of the supervision loop (default 30 days / 60 s / 600 s) |
| `wait_for_monte_carlo_output` | keep waiting for an external `3_monte_carlo/output` |
| `calculators`, `calculator_plugins` | per-calculator settings and extra interface modules |
| `<stage>` | settings of one stage (`"coll": {...}`), see [Stages](#stages-and-phases) |

### Stage instances

`run` may mix the same stage several times with different calculators:

```python
"run": [
    {"id": "coll_qe", "stage": "coll", "calculator": "espresso",
     "params": QE_PARAMS},
    {"id": "coll_vasp", "stage": "coll", "calculator": "vasp",
     "params": VASP_PARAMS},
    {"id": "relax_vasp", "stage": "relax", "calculator": "vasp",
     "coll_calculator": "espresso", "on_missing": "skip"},
    {"id": "mae_vasp", "stage": "mae", "calculator": "vasp",
     "coll_calculator": "vasp", "host": "cluster1", "max_parallel": 8},
],
```

| instance key | meaning |
| --- | --- |
| `id` | name of the instance (`--instances` selects it); unique |
| `stage`, `calculator` | **required** |
| `coll_calculator` / `coll_instance` | whose `2_coll` results this instance reads |
| `depends_on` | explicit instance dependencies |
| `on_missing`, `on_unsupported` | `error` (default) or `skip` when prerequisite results / stage support are missing |
| `reuse` | skip the instance when its output already exists |
| `host` | host name (or an inline host dictionary) for all steps of the instance |
| `execution`, `max_parallel`, `max_attempts`, `header` | per-instance execution plan |
| `steps` | phases of this instance, overriding the defaults |
| `work_dir` | explicit working directory of the instance |

Without the `id`/`stage` form, `run: ["coll", "mae"]` creates one instance per
stage, named after it.

## Running

| command | behaviour |
| --- | --- |
| `run_automag.py` | run the selected phases **once**, in order (default: `submit`) |
| `run_automag.py --wait` | run everything, polling the started jobs between steps; resumable |
| `run_automag.py --daemon` | same supervision loop, meant for a long-running process (nohup/systemd) |
| `run_automag.py --once` | advance until a step must wait for external programs, record the state, exit |
| `run_automag.py --status` | print the stored state |
| `run_automag.py --report` | detailed report: steps, hosts, running jobs and their ids |
| `run_automag.py --list-stages` | stages and phases |
| `run_automag.py --list-hosts` | the resolved hosts |
| `run_automag.py --dry-run` | build the plan and the stage inputs, execute nothing |

Selection and overrides: `--steps submit,postprocess`, `--stages coll,relax`,
`--instances coll_vasp,mae_vasp`, `--execution auto|inprocess|process|scheduler`,
`--state-file PATH`, `--poll-interval SECONDS`, `--wait-timeout SECONDS`,
`--fresh` (forget the stored state and start over).

### Supervised run (`--wait` / `--daemon`)

```
step 1  ──► launch its jobs ──► poll until they finish ──► verify results
                                  │
step 2  ◄─────────────────────────┘   (only when its dependencies are done)
```

Every step is a stage script (`1_submit.py`, `2_plot_results.py`, ...) started
by automag; its *calculations* are the individual VASP jobs.  The runner

1. starts a step when its dependencies are satisfied,
2. lets the step register the calculations it wants (a manifest in the stage
   directory) and starts them — up to `max_parallel` at a time,
3. polls their state (process id, or `squeue`/`qstat`), records every change in
   SQLite,
4. when all of them have results, continues with the next step.

With `--once` the same thing happens for exactly one pass: it starts what can
be started, records the job ids and exits; run the same command again later to
continue.  That makes it suitable for cron:

```cron
*/20 * * * * cd /work/project && /path/to/run_automag.py --once >> automag.log 2>&1
```

### State and resuming

The state file (`automag_state.sqlite`, `--state-file` to move it) holds one row
per **step** (`mae_vasp.submit`) and per **job**:

```
step      status, host, started_at, finished_at, config_hash, reset_at, detail
job       task_key, step, host, kind (process|slurm|pbs), external_id,
          launcher, workdir, command, items, attempts, status, detail
```

`--status` prints the summary, `--report` adds the running jobs and their
scheduler ids.  The database can be inspected with any SQLite client; it is the
only place where automag keeps its own bookkeeping.

Resuming is automatic: a step whose results exist is `done` and is skipped, a
job that is still alive is adopted instead of being started twice, a job that
disappeared without producing results is retried (`max_attempts`, default 3).

### Changing the settings in the middle of a run

Every step stores a hash of the settings it was started with.  Rerunning the
pipeline after editing `input_automag.py` compares the hashes:

* unchanged steps keep their status and are **not** recomputed;
* a step whose settings changed is reset: its recorded jobs and results are
  dropped, the affected calculations start again (a stale result on disk is not
  reused);
* a changed parameter that also changes the directory name (e.g. another
  `encut`) resets the step too, but only the calculations of the new parameter
  set are started — the old directories stay on disk and their results are not
  reused.

`--fresh` ignores the whole stored state and rebuilds the pipeline from
scratch.

## Where a step runs

Each stage instance picks a **host**; the host describes the machine and how
commands reach it:

```python
PIPELINE = {
    "hosts": {
        "local":   {"transport": "local"},
        "cluster": {"transport": "ssh",
                    "host": "login.cluster.example", "user": "me",
                    "identity_file": "~/.ssh/id_ed25519",
                    "automag_home": "/shared/automag",   # checkout on that machine
                    "python": "python3",
                    "scheduler": "slurm"},               # slurm | pbs | none
        "container": {"transport": "shell",              # commands through a wrapper
                      "wrapper": "docker exec -i automag"},
    },
    "run": [
        {"id": "coll",  "stage": "coll", "calculator": "vasp", "host": "cluster"},
        {"id": "mae",   "stage": "mae",  "calculator": "vasp", "host": "local"},
        {"id": "relax", "stage": "relax", "calculator": "vasp",
         "host": {"transport": "ssh", "host": "gpu.example", "user": "me",
                  "automag_home": "/opt/automag", "scheduler": "pbs"}},
    ],
}
```

| transport | how a command is executed |
| --- | --- |
| `local` (default) | `subprocess` on this machine |
| `ssh` | the same command prefixed with `ssh` (`BatchMode=yes`; keys must work) |
| `shell` | through `wrapper` on this machine, e.g. inside a container |

`workdir`, `automag_home`/`automag_exe`, `python`, `environment` (`module load`,
`source ...`), `ssh_options`, `port` are all optional; `--list-hosts` prints the
resolved table, and the resolved host of every step is recorded in the state
file.

### How a step is executed (smart default)

The *stage script* and its *calculations* are two different things:

* the stage script (`2_plot_results.py`, ...) runs either **in the automag
  process** (`inprocess`) or as a **separate process** (`process`).  `auto`
  (default) picks `process`, and falls back to `inprocess` on a machine that
  cannot afford another process.  A step on a remote host always runs as a
  separate process there.  `execution: "scheduler"` makes the script itself a
  batch job — never chosen automatically;
* the calculations are started with `sbatch`, `qsub` or as plain processes.
  `auto` (default) looks at the job header of the step: `#SBATCH` directives →
  SLURM, `#PBS` → PBS/Torque, otherwise a local process.  A host with
  `"scheduler": "slurm"|"pbs"` selects the scheduler for all its steps.

That is the "smart" behaviour the input file relies on: a header with batch
directives means batch jobs, no directives mean plain processes, and the host
decides *where*.

### Parallelism

Every calculation of a step is an independent job, so the calculations of a
step run in parallel up to the limit:

```python
"execution": {"max_parallel": 8},          # default for every step
"coll": {"max_parallel": 16},              # per stage
{"id": "mae_vasp", ..., "max_parallel": 4} # per instance
{"id": "mae_vasp", ..., "execution": {"mode": "process",
                                      "phases": {"submit": {"max_parallel": 2}}}}
```

Dependencies between calculations are respected automatically.  In `4_mae` the
reference calculation `z` finishes before the θ/φ scan starts, and the scan
finishes before the `MAE_curve` rotations — no manual ordering is needed.  A
calculation that dies without producing a result is retried (`max_attempts`,
default 3); after that the step continues with the results that do exist and
prints how many are missing, which is what you want when a single magnetic
configuration refuses to converge.  If the missing results are needed, the
analysis step fails with the list of failures and the path to its log.

Logs of the steps are in `<stage dir>/.automag_logs/`.

## Stages and phases

| stage | directory | phases |
| --- | --- | --- |
| `conv_tests` | `0_conv_tests` | `submit`, `postprocess` |
| `lin_response` | `1_lin_response` | `submit`, `postprocess` |
| `coll` | `2_coll` | `submit`, `postprocess` |
| `relax` | `2_relax` | `submit`, `postprocess` |
| `monte_carlo` | `3_monte_carlo` | `submit`, `run`, `espins`, `postprocess` |
| `mae` | `4_mae` | `submit`, `postprocess`, `axes`, `curve`, `curve_postprocess` |

`--wait` runs the phases in the right order, `2_coll` post-processing included,
so the dependent stages see the most stable configuration.  With `--steps` you
can run any subset by hand.

### `conv_tests`

Convergence of `ENCUT`, and of `SIGMA` + k-mesh:

```python
"conv_tests": {"mode": "kgrid",            # or "encut"
               "params": {...},
               "encut_values": range(500, 700, 50),
               "sigma_values": [0.05, 0.10],
               "kpts_values": range(10, 40, 10),
               "magnetic_atoms": ["Fe"],
               "configuration": "fm"}
```

`mode="encut"` sweeps `ENCUT`, `mode="kgrid"` sweeps `SIGMA` and `kpts`.
Defaults: ENCUT 500–1000 eV step 10, SIGMA 0.05–0.20 eV step 0.05, kpts 20–100
Å⁻¹ step 10.  `postprocess` plots energy vs parameter and prints the values
converged to 1 meV/atom w.r.t. the most accurate result.

### `lin_response`

Hubbard U by linear response.  Needs the dummy-atom POTCAR trick: the POTCAR of
the perturbed species is placed in the directory of `dummy_atom` of your
pseudopotential library.

```python
"lin_response": {"dummy_atom": "Zn", "dummy_position": 0,
                 "perturbations": [-0.1, -0.05, 0.0, 0.05, 0.1],
                 "params": {...}, "magnetic_atoms": ["Fe"]}
```

Produces `CalcFold/charges.txt` and the U plot; `U = 1/X - 1/X0`.

### `coll`

```python
"coll": {"supercell_size": 1,
         "spin_values": {"Fe": [4]},        # optional: HS, or [HS, LS]
         "magnetic_atoms": ["Fe"],          # optional: defaults to transition metals
         "lower_cutoff": 0,                 # optional: discard quenched moments
         "enumeration_mode": "enumlib",     # or "pymatgen"
         "params": {...}}
```

`enumlib` (default) uses the Wyckoff analysis driver; `pymatgen` uses
`MagneticStructureEnumerator` with `enumeration_strategies`
(`ferromagnetic`, `antiferromagnetic`, `antiferromagnetic_by_motif`,
`ferrimagnetic_by_motif`, `ferrimagnetic_by_species`, `nonmagnetic`) and an
optional `max_orderings`.  Both write the same `2_coll/trials_*/` layout
(`settingNNN.vasp` + `configurationsNNN.txt`), so later stages work with either.
`spin_values` and `magnetic_atoms` may be omitted — pymatgen's
`DEFAULT_MAGMOMS` are used.

`postprocess` writes the energy histogram, the table of configurations and
`..._most_stable_configuration.txt`, which the dependent stages read.

### `relax`

Re-relaxes the most stable configuration (`configuration: "lowest"`, or a
specific name), using the converged moments from `2_coll`:

```python
"relax": {"configuration": "lowest",
          "relax_params": {"ibrion": 2, "isif": 3, "nsw": 100, "ediffg": -0.01},
          "update_structure": False,      # True: copy the result into geometries/
          "params": {...}}
```

Result: `2_relax/relaxed_<state>.vasp` (and `.json`).

### `monte_carlo`

```python
"monte_carlo": {"configuration": "lowest",
                "cutoff_radius": 4.1,
                "control_group_size": 0.4,
                "append_coupling_constants": True,
                "magnetic_atoms": ["Fe"],
                "mc_backend": "vampire",          # or "espins"
                "mc_timeout": None}               # seconds
```

`submit` computes the exchange couplings (with the PCC of the Heisenberg model
as a quality measure) and writes the VAMPIRE unit cell / ESpinS inputs;
`run` executes the Monte Carlo code (`vampire_run_command`,
`espins_run_command`, `vampire_input_file`, `vampire_material_file`,
`espins_output_file`, `espins_output_converter`, `vampire_dir`, `espins_dir`,
...); `postprocess` fits the magnetization curve and prints Tc and the critical
exponent.  `wait_for_monte_carlo_output` keeps the old behaviour when you run
VAMPIRE by hand.

### `mae`

```python
"mae": {"configuration": "lowest",
        "Nph": 8, "Nth": 4, "N_MAE": 8,      # θ/φ scan and rotation curve
        "standardize_cell": True,
        "use_primitive_cell": True,
        "spin_values": {"Fe": [4]},
        "params": {...}}
```

Phases: `submit` (reference `z` + θ/φ scan), `axes` (easy/hard axis from the
scan), `curve` (rotations around the easy axis), `curve_postprocess` (K1, K2,
M0, (BH)max).  Output: `MAE_theta_phi_<conf>_U…_J…_K…_EN….txt`,
`MAE_curve_…txt`, plots and `.npy` arrays in `4_mae/outputs_<conf>/`.

#### Convergence study over the calculation parameters

Any sweepable parameter may be given as a **list** instead of a single value:

```python
"mae": {"Nph": 4, "Nth": 2, "N_MAE": 4,
        "params": {"encut": [400, 500, 600],   # 3 sets
                   "kpts": [10, 20],           # ×2 → 6 MAE calculations
                   "sigma": 0.1,               # fixed
                   "ldauu": [4.0, 0],          # NOT a sweep (per species)
                   ...}}
```

Sweepable keys: `encut`, `kpts`, `sigma`, `ediff`, `ediffg`, `prec`, `nelm`,
`ncore`, `lreal`, `ismear`, `kgamma`, `kspacing`, `lasph`, `amix`, `bmix`,
`nsw`, `ibrion`, `isif`.  Per-species lists (`ldauu`, `ldauj`, `ldaul`,
`setups`, `saxis`, `magmom`, ...) and an explicit `kpts` mesh (`[a, b, c]`) keep
their meaning.

Every combination becomes an independent MAE calculation in its own directory
(`..._K20_EN400_S01`, `..._K10_EN500_S02`, ...) with its own reference `z`, and
they are scheduled in parallel like any other calculations.  The analysis
writes one table per set plus

```
4_mae/MAE_convergence_<conf>.txt        MAE per set + relative change
4_mae/MAE_curve_convergence_<conf>.txt  MAE and K1 per set + relative change
```

```
set   parameters  MAE (MJ/m^3, per cell)
----  ----------  ----------------------
_S01  encut=400   47.736179
_S02  encut=500   47.736179

relative change of MAE (MJ/m^3, per cell) between successive sets:
  encut=400 -> encut=500: 0.00 %
```

`sweep_sets.json` in the stage directory records which sets belong to the
current run, so results of an earlier parameter set are never mixed in.

## Calculators and interfaces

The pipeline never branches on the calculator name; it asks the registry in
`calculators/`:

* built in: `vasp` (reference), `espresso` (alias `qe`), `fhi-aims` (alias
  `aims`), the generic `ase` interface and `prophet_spin` (alias `prophet`,
  needs the external model - not used by the examples);
* any ASE calculator: `"calculator": "ase:ase.calculators.espresso.Espresso"`;
* a plugin: a module with a `calculators.base.CalculatorInterface` subclass
  decorated with `calculators.register`, listed in `calculator_plugins` (or in
  `AUTOMAG_CALCULATORS`).

An interface implements `write_input`, `run_command` and optionally
`read_energy`, `read_magnetic_moments`, `read_structure`, `is_converged`,
`result_writer`; capability flags (`supports_ldau`, `supports_noncollinear`,
`supports_relaxation`) turn unsupported modes into clear errors or a skip
(`on_unsupported: "skip"`).  See `docs/AGENT_GUIDE.md` for a template, and
`examples/qe`, `examples/fhi-aims` for complete projects.

## External codes

| code | used by | default command |
| --- | --- | --- |
| `enumlib` | `coll` | `enum.x`, `makeStr.py` |
| `VAMPIRE` | `monte_carlo` (`vampire`) | `vampire` |
| `ESpinS` | `monte_carlo` (`espins`) | `mc.x` |

```bash
./build.sh --external                  # download + compile into external/
./build.sh --external --vampire-path /opt/vampire/vampire-serial \
                      --espins-path  /opt/espins/mc.x   # register what you have
```

The result is `external/codes.json`, read at run time by `common/external.py`
(pre/post commands such as `module load` included).  `AUTOMAG_EXTERNAL_CODES`
points at another registry, `PIPELINE["external_codes"]` overrides single
entries.  A code that is missing only disables the step that needs it, with an
explanatory message.  Details: [docs/INSTALL.md §6](docs/INSTALL.md).

## The bundled Prophet models

`models/prophet-v1-mag.pt` (energy model) and `models/prophet-v1-mag-head.pt`
(moment predictor) live in the automag checkout and are copied into the frozen
bundle, so `prophet_spin` works out of the box - no absolute paths in the
input file:

```python
"params": {"model_path": None, "device": "cpu"}   # None -> the bundled model
```

`PROPHET_MODEL` / `PROPHET_HEAD` / `PROPHET_DEVICE` override the defaults, and
an explicit `model_path` in the stage params always wins.  Without
`predict_moments` only the energy model is used; with it the head model predicts
the moments first (see `examples/prophet_moments/`).

The released weights are isotropic in the moment direction (energies of
`M || z` and `M || x` differ by less than 1e-6 eV), so `4_mae` with this model
returns MAE = 0.  That is why the default plan written by the bootstrap and the
`prophet/` example run `2_coll` only; `examples/prophet_mae/` shows how to run
the stages anyway (template for a direction-sensitive model) and VASP gives the
physical anisotropy.

## Paths and environment

| variable | meaning |
| --- | --- |
| `AUTOMAG_CODE_PATH` | where the sources live (set automatically) |
| `AUTOMAG_WORKDIR` | data root when the input file has no `workdir` |
| `AUTOMAG_PATH` | legacy: `2_coll`/`CalcFold` data root, used by the per-stage workflow |
| `AUTOMAG_EXTERNAL_CODES` | `codes.json` or the directory containing it |
| `AUTOMAG_CALCULATORS` | extra calculator plugin modules |
| `VASP_PP_PATH` | VASP pseudopotential library (for VASP runs) |
| `AUTOMAG_EXE` | compiled binary, used by `examples/*/run_pipeline.sh` |

## Build the executable

The Python sources need no compilation - this is only for a portable bundle.
Prerequisites: `pyinstaller` (and `ninja`, which `build.sh` uses):

```bash
conda install -c conda-forge ninja -y && python -m pip install pyinstaller
```

```bash
./build.sh                     # single file   -> dist/run_automag
./build.sh --onedir            # directory     -> dist/run_automag/ (recommended)
./build.sh --all               # external codes + binary
./build.sh --external          # only enumlib / VAMPIRE / ESpinS -> external/
python tools/build_executable.py --dry-run       # show the PyInstaller command
python tools/build_executable.py --no-prophet    # leave torch/prophet_spin out
python tools/build_executable.py --no-models     # prophet inside, checkpoints outside
```

Two ready-made variants differ only in where the Prophet checkpoints come from:

| build | size | prophet |
| --- | --- | --- |
| `--onedir --name run_automag` | ~5.7 GB | inside: `prophet_spin`, torch and the checkpoints (`_internal/models/`), nothing to configure |
| `--onedir --name run_automag-extmodels --no-models` | ~5.4 GB | inside, but the checkpoints stay outside: `PROPHET_MODEL` / `PROPHET_HEAD` must point at the `.pt` files |
| `--onedir --name run_automag-noprophet --no-prophet --no-models` | ~500 MB | absent: VASP / Quantum ESPRESSO / FHI-aims / ASE only |

The first two understand `calculator: prophet_spin`.  The `--no-models` variant
keeps the 275 MB of weights out of the bundle (handy when one copy of `models/`
is shared over NFS); the `--no-prophet` variant is the small DFT-only build -
it starts in a directory with just a structure too, but says plainly that it
cannot run the generated prophet plan.

What ends up inside:

| part | notes |
| --- | --- |
| Python + libraries | its own interpreter, no conda needed at run time |
| stage scripts | `0_conv_tests/`, `2_coll/`, `2_relax/`, `3_monte_carlo/`, `4_mae/`, ... (without the example data) |
| pymatgen / ase data | `periodic_table.json`, `ase.io.db`, ... |
| `models/*.pt` | both Prophet checkpoints (275 MB) |
| `prophet_spin` + torch | bundled by default; `--no-prophet` removes them |
| `external/codes.json` | the registry of enumlib / VAMPIRE / ESpinS, when present |

Sizes: ~490 MB without `prophet_spin`, ~5.7 GB with it (torch).

## Run the compiled version

The bundle is self-contained: it brings its own Python and libraries and runs
**without the conda environment it was built in**, and without the sources:

```bash
# on any machine with the same OS/arch, no python and no conda needed
env -i PATH=/usr/bin:/bin ./dist/run_automag/run_automag --dry-run --wait
env -i PATH=/usr/bin:/bin ./dist/run_automag/run_automag --wait
```

In a project directory the launchers take care of that:

```bash
cd my_project
AUTOMAG_EXE=/path/to/dist/run_automag/run_automag ./run_pipeline_compiled.sh
AUTOMAG_EXE=/path/to/dist/run_automag/run_automag ./run_pipeline_slurm_compiled.sh
```

The Prophet models travel with the bundle (`<bundle>/_internal/models/`), so
`calculator: prophet_spin` works with no extra file and no `model_path` - keep
`params` empty or set `"device": "cuda"` when the build has CUDA.  Building with
`--onedir` is recommended: a single-file build unpacks ~6 GB to a temporary
directory at every start.

See [docs/INSTALL.md §5](docs/INSTALL.md) for the full build notes.

## Examples and tests

```bash
cd examples/vasp && python ../../run_automag.py --dry-run --wait
```

`examples/` contains `vasp`, `qe`, `fhi-aims`, `pymatgen_enumeration`,
`prophet` (Prophet-V1-Mag instead of DFT), the mixed setups
`mixed_fhi_aims_vasp` / `mixed_qe_vasp`, and `remote_vasp`,
`remote_fhi_aims`, `remote_qe`, which run `2_coll` on another machine
(`PIPELINE["hosts"]`); `examples/README.md` explains them and their launchers.

DFT-free test suites (`tools/`):

| test | what it checks |
| --- | --- |
| `smoke_test.py` | calculators, `reuse`/`skip` policies, external-code wiring |
| `execution_test.py` | `--once` barrier, resuming, parallel limit, config change reset, SLURM/PBS launcher, remote host |
| `mae_chain_test.py` | the complete `1_submit → 2_plot → 3_submit → 4_plot` chain and the parameter sweep, with synthetic VASP output |

## Legacy per-stage workflow

The old way of running each stage by hand still works: edit `input.py` in
`0_conv_tests/`, `2_coll/`, ... and run `1_submit.py` / `2_plot_results.py`
with `AUTOMAG_PATH` and (optionally) FireWorks pointing at your MongoDB
launchpad.  The stage scripts are the same ones the pipeline starts, so results
are interchangeable.  New projects should use the pipeline — it replaces the
manual `qlaunch`/`sbatch` bookkeeping described in older revisions of this
README.

## License and citation

See the repository for the license.  If you use automag in a publication,
please cite the project as described there.