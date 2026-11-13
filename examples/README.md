# Automag examples

Self-contained examples for running the pipeline in the automatic mode.  Each
directory contains:

- `structure.vasp` - the input structure (Fe4O6, alpha phase);
- `input_automag.py` - the controlling file (stages, instances, settings).
  `input_automag.yaml` / `.json` work too (priority `.py` → `.yaml` → `.json`);
  the launcher picks whichever exists;
- the launchers, one per way of running the pipeline (see below).

The directories are meant to be copied anywhere.  Every launcher changes into
its own directory first, so all results (`0_conv/`, `2_coll/`, `2_relax/`,
`3_monte_carlo/`, `4_mae/`, `CalcFold/`, `geometries/`, `automag_state.sqlite`)
are written next to `input_automag.py` and the automag checkout is never
modified.

## Launchers

Six scripts, two independent choices: **how** the pipeline is run (one blocking
command, one pass at a time, or a SLURM job) and **what** runs it (the Python
sources or the compiled binary).

| script | runs the supervisor | where | returns |
| --- | --- | --- | --- |
| `run_pipeline.sh` | Python sources (`$AUTOMAG_HOME/run_automag.py`) | current shell | when the whole pipeline is finished |
| `run_pipeline_compiled.sh` | compiled binary (`$AUTOMAG_EXE`) | current shell | when the whole pipeline is finished |
| `run_pipeline_once.sh` | sources or binary | current shell | after **one pass**, work may remain |
| `run_pipeline_once_loop.sh` | sources or binary | current shell | when finished, looping over passes |
| `run_pipeline_slurm.sh` | Python sources | SLURM job (self-submitting) | job id, immediately |
| `run_pipeline_slurm_compiled.sh` | compiled binary | SLURM job (self-submitting) | job id, immediately |

All of them accept the same pipeline options and print the chosen binary and
controlling file.

`AUTOMAG_HOME` (a checkout) and `AUTOMAG_EXE` (a binary) are alternatives: the
two `*_compiled.sh` scripts require the binary and never fall back to the
sources, so a typo cannot silently start a different installation.

The compiled build contains VASP, Quantum ESPRESSO, FHI-aims and the generic
ASE calculators; `prophet_spin` is not part of it by default (it would add
several GB of torch/CUDA libraries) and is not used by these examples.

## The two launch styles

### One command that waits

```bash
cp -r /path/to/automag/examples/qe ~/my-qe-run
cd ~/my-qe-run
export AUTOMAG_HOME=/path/to/automag      # source checkout
# or: export AUTOMAG_EXE=/shared/bin/run_automag   # compiled binary
export AUTOMAG_CONDA_ENV=automag
./run_pipeline.sh
```

`--wait` submits each batch of calculations, polls the jobs until they finish,
runs the analysis, submits the next batch and so on.  Nothing has to be
restarted by hand; the process keeps the state in `automag_state.sqlite`, so if
it is killed (logout, wall time, Ctrl-C) the same command continues where it
stopped.

### Manual passes: one pass per command, or a loop

The same pipeline can be advanced one pass at a time.  A pass starts whatever
can be started, records the job ids and exits:

```bash
./run_pipeline_once.sh          # pass 1: submits the first calculations
# ... the cluster runs them ...
./run_pipeline_once.sh          # pass 2: analysis + the next calculations
./run_pipeline_once.sh          # ... until "every step is complete"
```

`run_pipeline_once_loop.sh` automates exactly that: one `--once` pass, then a
sleep, then the next pass, until the pipeline reports that it is complete.

```bash
export AUTOMAG_ONCE_INTERVAL=300          # seconds between passes (default 60)
./run_pipeline_once_loop.sh

# a bounded number of passes (useful to watch it work)
AUTOMAG_ONCE_MAX_ROUNDS=3 ./run_pipeline_once_loop.sh
```

Because every pass only looks at the state file, this is also what a cron job
does - no long-running process is needed at all:

```cron
*/20 * * * * cd /work/project && /path/to/automag/examples/qe/run_pipeline_once.sh >> automag.log 2>&1
```

`run_pipeline.sh` and `run_pipeline_once_loop.sh` do the same work; the first
sleeps inside one process, the second between processes.  Use the first when a
single blocking command is convenient, the second when the supervisor must not
stay in one process (short wall-time limits, reboots, manual control).

### Running the supervisor under SLURM

The `*_slurm*.sh` scripts turn *the supervisor itself* into a batch job, so it
survives the logout and does not sit on the login node.  The calculations it
launches are ordinary SLURM jobs - automag takes the scheduler from the
`#SBATCH` directives in the step header (or `"scheduler": "slurm"` in
`PIPELINE["hosts"]`).

```bash
export AUTOMAG_HOME=/path/to/automag
./run_pipeline_slurm.sh                # prints the job id, returns immediately
squeue -u $USER                        # find the driver
tail -f automag-driver-*.out           # follow it

# resources are in the #SBATCH header; override them without editing the file
AUTOMAG_SLURM_OPTS="--partition=compute --account=myproj --time=2-00:00:00" \
    ./run_pipeline_slurm.sh

# pipeline options are forwarded to the job
./run_pipeline_slurm.sh --dry-run
./run_pipeline_slurm_compiled.sh --stages coll --steps submit
```

The driver job must be allowed to call `sbatch` itself (it submits the
calculations).  If your site forbids that, use `run_pipeline_once_loop.sh` from
cron instead - the state file makes the two approaches interchangeable.

### Environment variables understood by the launchers

| variable | used by | meaning |
| --- | --- | --- |
| `AUTOMAG_HOME` | `run_pipeline.sh`, `run_pipeline_slurm.sh`, `*_once*.sh` | automag checkout (default: this example's `../../`) |
| `AUTOMAG_EXE` | the compiled launchers | compiled `run_automag` binary (a binary on `PATH` is used automatically) |
| `AUTOMAG_CONDA_ENV` | the source launchers | conda env with the dependencies (default `automag`; empty string = do not activate anything) |
| `AUTOMAG_CONDA_SH` | the source launchers | `conda.sh` to source when `conda` is not on `PATH` (typical inside a batch job) |
| `AUTOMAG_ONCE_INTERVAL` | `run_pipeline_once_loop.sh` | seconds between passes (default 60) |
| `AUTOMAG_ONCE_MAX_ROUNDS` | `run_pipeline_once_loop.sh` | stop after N passes (default 0 = until finished) |
| `AUTOMAG_PRE_COMMAND` | the SLURM launchers | shell line executed at the start of the job, e.g. `"module load python/3.12"` |
| `AUTOMAG_SLURM_OPTS` | the SLURM launchers | extra `sbatch` options (`--partition`, `--account`, `--time`, ...) |
| `AUTOMAG_EXTERNAL_CODES` | all | registry of enumlib/VAMPIRE/ESpinS; a `codes.json` next to the binary is picked up automatically |

Conda activation is done with `set +u` around it: conda's own deactivation
hooks are not `set -u` safe, and a launcher that dies during `conda activate`
is impossible to debug from a batch log.

## Before the first run

Edit the placeholders in `input_automag.py`:

- the SLURM header (`--partition`, `--ntasks`, `--gres`, ...);
- `calculator_command` (module names, executable, MPI launcher);
- calculator-specific paths (pseudopotentials for QE, `species_dir` for
  FHI-aims);
- `environment_activate` / `environment_deactivate` used inside the job
  scripts.

Extra arguments are passed through to the pipeline:

```bash
./run_pipeline.sh --stages coll
./run_pipeline.sh --stages monte_carlo --steps submit
./run_pipeline.sh --instances coll_qe
./run_pipeline.sh --dry-run            # print the plan only
./run_pipeline.sh --wait-timeout 604800
```

## What each example runs

| calculator(s) | directory | stages |
| --- | --- | --- |
| VASP | `vasp/` | `conv_tests`, `coll`, `relax`, `monte_carlo`, `mae` |
| Quantum ESPRESSO | `qe/` | `coll`, `monte_carlo` |
| FHI-aims | `fhi-aims/` | `coll`, `monte_carlo` |
| VASP + pymatgen enumeration | `pymatgen_enumeration/` | `coll` (`MagneticStructureEnumerator`), `relax`, `monte_carlo` |
| FHI-aims → VASP | `mixed_fhi_aims_vasp/` | `coll(aims)`, `mc(aims)`, `mae(vasp, coll=aims)` |
| QE → VASP | `mixed_qe_vasp/` | `coll(qe)`, `mc(qe)`, `mc(vasp, coll=qe)`, `mae(vasp, coll=qe)` |
| Prophet-V1-Mag (no DFT) | `prophet/` | `coll` (`monte_carlo` degenerates at this cell size) |
| Prophet, moments predicted first | `prophet_moments/` | `coll` (two-step: moments, then energy) |
| Prophet for the anisotropy | `prophet_mae/` | `coll`, `mae` (MAE ~ 0 - not in the default plan, see below) |
| Prophet (local) + VASP (cluster), original and relaxed cell | `prophet_then_vasp/` | 13 instances: `coll`/`monte_carlo`/`mae` twice per calculator, `relax` with VASP |
| VASP on a remote cluster | `remote_vasp/` | `coll` (on the host `cluster`) |
| FHI-aims on a remote cluster | `remote_fhi_aims/` | `coll` (on the host `cluster`) |
| QE on a remote cluster | `remote_qe/` | `coll` (on the host `cluster`) |

`coll` computes the energies of the generated magnetic configurations and
finds the ground state.  By default it uses automag's own enumlib driver;
`pymatgen_enumeration/` instead sets `enumeration_mode = "pymatgen"` to use
`pymatgen`'s `MagneticStructureEnumerator` (strategies such as
`antiferromagnetic_by_motif`, `ferrimagnetic_by_species`, `max_orderings`).
Neither needs `spin_values`/`magnetic_atoms`: missing moments come from
pymatgen's `DEFAULT_MAGMOMS`.  Both write the same `trials_*/` layout.

`monte_carlo` computes the Heisenberg coupling constants and writes the
VAMPIRE/ESpinS input; the `run` phase launches the Monte Carlo backend
automatically.  The Monte Carlo examples use `configuration = "lowest"`, which
resolves to the most stable configuration found by `coll`.

The VASP example additionally runs convergence tests, a structural
re-relaxation of the most stable configuration, and the magnetocrystalline
anisotropy.  `lin_response` is not enabled by default because it requires the
dummy-atom POTCAR trick described in the main README; add it to `run` when
needed.

`lin_response`, `relax`, and `mae` currently use VASP-specific features
(LDA+U perturbations, VASP relaxation tags, non-collinear SAXIS).  The
Quantum ESPRESSO and FHI-aims examples therefore run the calculator-independent
part of the workflow.  The interfaces can be extended without touching the
pipeline - see `calculators/` in the main README.

The `mixed_*` examples show how to combine calculators in one run.  Each
instance is a dict with `id`, `stage`, `calculator` and optionally
`coll_calculator` (whose `2_coll` results it reads), `reuse`, `on_missing` and
`on_unsupported`.  The `monte_carlo`/`mae` instances are skipped when their
input data is missing or their output already exists - the plan is printed by
`./run_pipeline.sh --dry-run`.

`prophet/` runs the magnetic search with the **Prophet-V1-Mag** model instead
of a DFT code: `prophet_spin` is an ASE calculator, so the whole `coll` stage
works without VASP/QE/FHI-aims.  Its `run` list contains only `coll` on
purpose - see "the Heisenberg model needs a solvable system" below.  It needs
`torch` + `prophet_spin` in the environment; the checkpoints ship with automag
in `models/` (`prophet-v1-mag.pt` + `prophet-v1-mag-head.pt`) and are bundled
into the compiled binary, so no absolute paths are needed.  The
`PROPHET_MODEL` / `PROPHET_HEAD` / `PROPHET_DEVICE` environment variables
override them (`cpu` by default).  Inference on one core is slow, so give the
job several cores (`--cpus-per-task`) or add `--gres=gpu:1` and use `cuda`.

A directory that contains only a structure file needs no input file: automag
writes a default `input_automag.py` - the plan of `prophet/` (`coll` with the
bundled model; `4_mae` is left out because the model is isotropic, see
`prophet_mae/` to add it deliberately) - and runs it.  Delete the generated
file to get it again.

`relax`, `mae` and `lin_response` are VASP-specific: with a spin-lattice model
the pipeline skips such instances with a warning (`on_unsupported`).

### Prophet in two steps: moments, then energy

`prophet_moments/` uses the two checkpoints of Prophet-V1-Mag the way they are
meant to be combined.  One calculation is

1. `prophet-v1-mag-head.pt` (`MagmomPredictor`) predicts the magnetic moment of
   every atom from the initial spin configuration written by `2_coll`;
2. `prophet-v1-mag.pt` (`SpinCalculator`) evaluates the energy for those
   predicted moments.

Enable it with two keys in the stage `params`:

```python
PROPHET_PARAMS = {
    # both checkpoints ship with automag in models/ - no absolute paths
    "model_path": None,
    "head_path": None,
    "device": "cpu",
    "predict_moments": True,
}
```

`results.json` of such a calculation holds `energy`, `magmoms` (as before) and
`magmoms_predicted`, so the post-processing sees the predicted moments as the
converged ones.  For Fe4O6 `predict_moments = True` gives Fe moments of
+/-4.3675 and changes the energy of the ground state from -66.748 to
-67.008 eV.  Without `predict_moments` the model is evaluated for the initial
moments, which is what `prophet/` does.

### Prophet and 4_mae

`prophet_mae/` runs `4_mae` on prophet energies: the interface carries the
rotated moments as three numbers per atom, so `4_mae` writes its files and
`MAE(theta, phi)` is produced as usual.  The released Prophet-V1-Mag weights
are, however, isotropic in the moment direction (energies of `M || z` and
`M || x` differ by less than 1e-6 eV), so the MAE is numerically zero and
`MAE_curve` is flat.  The example is the template for a direction-sensitive
model; with VASP the same stages give the physical anisotropy.

### Two instances of one stage: `structure` and `struct_suffix`

Every instance of `run` may override two settings that are otherwise global:

| key | meaning |
| --- | --- |
| `structure` | input geometry, a file in `<workdir>/geometries/` (default: the global `structure`) |
| `struct_suffix` | suffix of the result directories, `CalcFold/<formula><suffix>/<calculator>/` |

That is how `prophet_then_vasp/` runs `coll`, `monte_carlo` and `mae` twice per
calculator - once for the input cell and once for the relaxed one - without
the two runs overwriting each other.  A `relax` instance with
`"update_structure": True` writes its result to
`geometries/<structure stem>_relaxed.vasp`, which is exactly the file the second
pass reads; give the second `coll` a `depends_on` (or a `coll_instance`) so it
starts after the relaxation.  Because each calculator then has two coll
instances, every `monte_carlo`/`mae` must name its `coll_instance`.

### Example layout and file formats

Every example is a self-contained directory with the launchers, the structure
and the settings.  The new examples also come in three formats, one directory
per format - `py/`, `yaml/`, `json/` - holding the same settings:

```
examples/prophet_moments/
    py/input_automag.py
    yaml/input_automag.yaml
    json/input_automag.json
```

The loader tries `input_automag.py`, then `.yaml`, then `.json`, so each
variant directory runs with `./run_pipeline.sh` exactly like the others.
All `input_automag.py` files follow one layout: parameter dictionaries at the
top (`VASP_PARAMS`, `PROPHET_PARAMS`, `..._SUBMIT`), then a `PIPELINE` with an
explicit `run` list of instances in the stage order of `examples/vasp`
(`conv_tests`, `coll`, `relax`, `monte_carlo`, `mae`).

## Running a stage on another machine

The `remote_*` examples run **`2_coll` on a different machine** while the
pipeline itself stays where it was started.  The instance says where its work
belongs:

```python
"hosts": {
    "local":   {"transport": "local"},
    "cluster": {"transport": "ssh",
                "host": "login.cluster.example", "user": "myuser",
                "workdir": "/shared/myproject",     # project dir there
                "automag_home": "/shared/automag",  # automag checkout there
                "python": "python3",
                "scheduler": "slurm"},
},
"run": [
    {"id": "coll_vasp_remote", "stage": "coll", "calculator": "vasp",
     "host": "cluster",          # <-- this stage runs on the cluster
     "supercell_size": 1, "spin_values": {"Fe": [4]}, "params": VASP_PARAMS},
],
```

* `remote_vasp/`, `remote_fhi_aims/`, `remote_qe/` differ only in the
  calculator and the parameters; `python ../../run_automag.py --list-hosts`
  prints what was resolved (`cluster: ssh -> myuser@login.cluster.example
  (scheduler: slurm)`).
* The job header of the stage carries `#SBATCH` directives, so the
  calculations are submitted with `sbatch` *on the cluster*; with
  `"scheduler": "pbs"` automag uses `qsub` instead, and without a scheduler it
  runs them as plain processes there.
* Requirements: the project directory must be visible on both sides at the same
  path (shared filesystem - results are read back from it), an automag checkout
  or compiled binary must exist there, and `ssh` must work without a password
  (`BatchMode=yes`; `identity_file` may point at a key).
* `environment` (a list of shell lines such as `module load vasp/6.4.3`) and
  `automag_exe` cover the module system and a compiled installation on the
  remote side.
* The pipeline can mix hosts freely: one stage locally, the next on the
  cluster, a third on a second cluster.
* No second machine at hand?  Switch the host to the wrapper transport and the
  same code path runs locally:
  `"cluster": {"transport": "shell", "wrapper": "sh -c", "scheduler": "none"}`.
  `tools/remote_examples_test.sh` does exactly that and checks that the stage
  really ran on the configured host.

### How many exchange couplings the fit can support

`monte_carlo` fits the exchange couplings from the `coll` energies, so it needs
at least as many independent magnetic configurations as couplings.  A 10-atom
cell with 4 magnetic atoms (Fe4O6, `supercell_size = 1`) gives 5
configurations - 2 independent equations for ~6 unknowns at the default
`cutoff_radius = 4.1`.

The stage does not abort in that case: it refits with `cutoff_radius` reduced
by 0.1 A until the system is solvable, never below the nearest-neighbour shell:

```
WARNING: system of 2 independent equation(s) in 6 unknowns at cutoff_radius = 4.10 A.
         reducing cutoff_radius to 4.00 A and refitting
...
WARNING: cutoff_radius was reduced from 4.10 to 2.90 A so that 1 coupling(s)
         can be fitted with 2 independent equation(s).
```

The reduced value is written back into the stage `input.py`, so
`2_write_vampire_ucf.py` and `5_run_mc.py` use the same shell list.  A fit like
this one is degenerate - one coupling, a perfect-looking correlation - so treat
it as a placeholder: for a real exchange model use `"supercell_size": 2` (or 3)
in `coll`, which the `vasp` example does.  Only if even a single shell cannot
be fitted does the stage stop with `SYSTEM OF N EQUATION(S) IN M UNKNOWNS`.

## Calculator notes

### VASP

Set `calculator_command` to your VASP launcher (e.g.
`module load vasp/6.4.3; mpirun vasp_std`).  The example uses LDA+U for Fe and
`isym = 0` so that the generated magnetic configurations are not symmetrized
away.

### Quantum ESPRESSO

ASE writes `espresso.pwi`; run `pw.x -in espresso.pwi > espresso.pwo`.
Put the UPF files in `pseudo_dir` and keep the `pseudopotentials` dictionary
consistent with them.

Quantum ESPRESSO defines `starting_magnetization` per atomic species.  The
pipeline writes the generated moments to the ASE `Atoms` object, but
antiferromagnetic order within one species is only fully representable if the
structure distinguishes the sublattices (for example by using different
species labels).  For a structure with a single Fe species, the AFM
configurations may collapse to the same QE state; the FM/NM configurations are
always distinct.

### FHI-aims

ASE writes `geometry.in` and `control.in`; `geometry.in` contains the
`initial_moment` of every atom, so all generated magnetic configurations are
reproduced exactly.  Set `species_dir` to your FHI-aims species defaults and
`calculator_command` to `aims.x` (or an MPI launcher).

## Monte Carlo backend

The examples use VAMPIRE (`mc_backend = "vampire"`).  The `run` phase copies
the VAMPIRE input and material file from `3_monte_carlo/vampire_input` (or the
paths you set) and runs the configured command.

Both VAMPIRE and ESpinS can be downloaded and compiled by the build:

```bash
./build.sh --external          # enumlib + VAMPIRE + ESpinS
```

If a code is already installed, register it instead (nothing is downloaded):

```bash
./build.sh --external \
    --vampire-path /opt/vampire/vampire-serial \
    --vampire-pre-command "module load intel" \
    --vampire-post-command "module purge"
```

without this the default commands are `vampire` / `mc.x` on `PATH`.  A code
that fails to compile is reported with a warning and its step stops with the
same instruction; see `README.md` → “External codes”.

For ESpinS set:

```python
"mc_backend": "espins",
# "espins_run_command": "/path/to/mc.x",   # default: external/codes.json
# "espins_output_file": "result.dat",
# "espins_output_converter": "python convert_espins.py",
```

If the Monte Carlo code is not available, remove `"monte_carlo"` from `run` or
run only the earlier stages.

## Checking the examples without DFT

Every example plans without a DFT code installed:

```bash
cd examples/vasp && python ../../run_automag.py --dry-run --wait
cd examples/qe   && python ../../run_automag.py --dry-run --once
```

`--dry-run` writes the stage inputs but starts nothing, so it is a safe way to
check a launcher and the controlling file after editing them.
