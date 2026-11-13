# Installing and building automag

This document covers the three supported ways of getting a working automag
installation, plus the optional compilation steps.  Nothing here is
distribution-specific unless it is marked as such.

Automag consists of three independent parts - install only what you need:

| part | needed for | how |
| --- | --- | --- |
| **Python runtime** | running the pipeline (`run_automag.py`) | conda, venv or system packages |
| **external codes** (enumlib, VAMPIRE, ESpinS) | the `coll`, `monte_carlo` stages | `./build.sh --external` or register existing binaries |
| **single-file build** (ninja + PyInstaller) | `dist/run_automag` binary | `./build.sh` |

Everything except the optional single-file binary is plain Python: **automag
does not have to be compiled to run**.

---

## 0. Requirements at a glance

### Python

* **Python 3.10 or newer** (3.12 is what the pipeline is developed and tested
  with).  The code uses `X | None` type annotations at run time, so 3.9 and
  older fail immediately.
  *Attention:* the distribution Python of RHEL/Rocky/Alma 9 is **3.9** - see
  [§3](#3-system-packages-no-conda-no-venv).
* `pip`, `wheel`, `setuptools`.

### Python packages

The pipeline itself needs only these:

```
numpy  scipy  matplotlib  ase  pymatgen  spglib
```

Optional:

| package | needed for |
| --- | --- |
| `PyYAML` | `.yaml`/`.yml` input files (`.py` and `.json` work without it) |
| `FireWorks`, `pymongo` | the legacy FireWorks database workflow (`"use_fireworks": True`) |
| `pyinstaller`, `ninja` | building `dist/run_automag` |
| `prophet_spin` | optional (`requirements-prophet.txt`): the Prophet-V1-Mag calculator; the checkpoints ship in `models/` |

`requirements.txt` is the full historical list (it also contains `Flask`,
`paramiko`, `fabric`, `plotly`, ...).  Installing it is the simplest and safest
option; the pipeline does not use the extra entries.  `ssh` access to remote
hosts goes through the system `ssh` client, not through paramiko.

Versions known to work together (from a reference installation):

```
Python 3.12.0   numpy 2.2.6   scipy 1.14.1   pymatgen 2024.10.22
ase 3.29.0      matplotlib 3.9.2   spglib 2.5.0   PyYAML 6.0.2
FireWorks 2.0.4 ruamel.yaml 0.18.6
```

### External programs

| program | needed for | note |
| --- | --- | --- |
| `gfortran`, `make`, `git` | compiling **enumlib** | system compilers |
| `g++`, `make`, `git` | compiling **VAMPIRE** | ~40 min single core |
| `ifort`/`ifx` + MKL, `make`, `git` | compiling **ESpinS** | its `make.sys` is Intel-only |
| `sbatch` (SLURM) or `qsub` (PBS/Torque) | submitting batch jobs | only if you use a scheduler |
| `ssh` | running steps on a remote host | password-less keys recommended |

### VASP (the reference calculator)

* `VASP_PP_PATH` pointing at the pseudopotential library, with `potpaw_PBE`
  (and `potpaw_LDA` if you use it):
  `export VASP_PP_PATH=/path/to/pp`.
* A VASP command that works on the compute nodes, e.g.
  `"calculator_command": "module load vasp/6.4.3; mpirun vasp_std"` in the
  pipeline input file.  Automag writes the inputs itself and only runs that
  command.

---

## 1. conda / mamba (recommended)

Conda gives you a modern Python and the compilers without touching the system
installation, and it is the only option that works out of the box on RHEL 9
with its Python 3.9.

```bash
# 1. create the environment (any name, "automag" below)
conda create -n automag python=3.12 -y

# 2. activate it
source /path/to/miniconda3/bin/activate automag
#   or:  conda activate automag   (after `conda init bash`)

# 3. get the sources
git clone https://github.com/Kurufinve/automag.git
cd automag

# 4. Python dependencies
python -m pip install --upgrade pip wheel setuptools
python -m pip install -r requirements.txt

# 5. only for the single-file build
conda install -c conda-forge ninja -y
python -m pip install pyinstaller
```

Verify:

```bash
python -c "import numpy, scipy, matplotlib, ase, pymatgen, spglib; print('ok')"
python run_automag.py --list-stages
```

Notes:

* `pymatgen` and `spglib` ship binary wheels for common platforms; if pip
  starts compiling them you are on an unsupported Python/platform combination -
  use a different Python version instead of fighting the build.
* If you also need to *compile* the external codes inside the environment
  (`conda install -c conda-forge gfortran gxx make git cmake`), be aware that
  `tools/build_external.py` deliberately removes the conda toolchain from
  `PATH` and puts `/usr/bin` first: a conda `g++` linked against a newer glibc
  produces binaries that do not run on the login node.  Therefore the system
  compilers (`gcc`, `g++`, `gfortran`) must also be installed - see §3.
* Environment variables (optional, but convenient).  Put these in `~/.bashrc`:

```bash
export VASP_PP_PATH=/path/to/pp
export AUTOMAG_EXTERNAL_CODES=/path/to/automag/external/codes.json   # optional
```

`AUTOMAG_PATH` is only needed for the legacy per-stage workflow (`cd 2_coll &&
python 1_submit.py`); the unified pipeline knows its own location.

---

## 2. venv without conda

Use this when the machine already has Python 3.10+ and you do not want conda.
Nothing else is required: the pipeline is pure Python, and the scientific
packages are available as wheels.

```bash
# check the interpreter first - it must be 3.10 or newer
python3 --version

git clone https://github.com/Kurufinve/automag.git
cd automag

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip wheel setuptools
python -m pip install -r requirements.txt

# only for the single-file build
python -m pip install pyinstaller
```

`ninja` is *not* a Python package you need to pip-install for automag's
single-file build; either install the `ninja` binary from your distribution
(`dnf install ninja-build` / `apt install ninja-build`) or create the build
without it:

```bash
python tools/build_executable.py            # PyInstaller directly, no ninja
```

Verify:

```bash
python -c "import numpy, scipy, matplotlib, ase, pymatgen, spglib; print('ok')"
python run_automag.py --list-stages
python tools/smoke_test.py                  # fast end-to-end test, no DFT
```

To activate the environment later: `source /path/to/automag/.venv/bin/activate`.

If `python3 -m venv` fails with `ensurepip is not available`, install the venv
support of your distribution: `apt install python3-venv` on Debian/Ubuntu,
`dnf install python3-pip` on RHEL 9 (the `venv` module itself is part of the
system `python3` package there).

---

## 3. System packages (no conda, no venv)

Everything is installed into the system Python.  It works, but it is the least
isolated option: a conflicting `numpy`/`pymatgen` version breaks automag and
everything else on the machine.  Prefer §1 or §2 unless the machine is
dedicated.

### 3.1 RHEL 9 / Rocky 9 / Alma 9

The system Python is **3.9**, which is too old.  Install a newer interpreter
from the distribution (`python3.11`/`python3.12` are in AppStream on current
RHEL 9 point releases; check with `dnf list available 'python3.1*'`):

```bash
sudo dnf install -y python3.12 python3.12-pip python3.12-devel
sudo dnf install -y gcc gcc-c++ gfortran make git cmake
sudo dnf install -y ninja-build          # only for ./build.sh

git clone https://github.com/Kurufinve/automag.git
cd automag
python3.12 -m pip install --upgrade pip wheel setuptools
python3.12 -m pip install -r requirements.txt
```

If no newer interpreter is packaged, use conda (§1) or build Python from
source; do not try to run automag on 3.9.

### 3.2 Debian / Ubuntu

```bash
sudo apt update
sudo apt install -y python3 python3-venv python3-pip python3-dev \
                    build-essential gfortran git cmake ninja-build

git clone https://github.com/Kurufinve/automag.git
cd automag
python3 -m pip install --user --upgrade pip wheel setuptools
python3 -m pip install --user -r requirements.txt
```

`--user` installs into `~/.local`; make sure `~/.local/bin` is on `PATH`.

### 3.3 Any distribution, without root

```bash
python3 -m pip install --user -r requirements.txt
```

Same interpreter-version requirement (≥ 3.10).

---

## 4. Verifying the installation

Run these from the automag checkout.  All of them work without any DFT code.

```bash
# 1. everything compiles
python -m py_compile run_automag.py pipeline/*.py calculators/*.py \
      common/*.py 0_conv_tests/*.py 1_lin_response/*.py 2_coll/*.py \
      2_relax/*.py 3_monte_carlo/*.py 4_mae/*.py

# 2. the plan of every example is built (no command is executed)
for d in examples/*/; do
    (cd "$d" && python ../../run_automag.py --dry-run --wait >/dev/null) || echo "FAILED $d"
done

# 3. unit-ish end-to-end tests (~1 min each, no DFT)
python tools/smoke_test.py          # calculators, reuse/skip, external codes
python tools/execution_test.py      # --once, resume, SLURM/PBS, remote host
python tools/mae_chain_test.py      # full 1-4 stage chain + parameter sweep

# 4. what the pipeline sees
python run_automag.py --list-stages
cd examples/vasp && python ../../run_automag.py --list-hosts
```

---

## 5. Building the executable

Optional.  It bundles Python, the dependencies, the stage scripts, the
`geometries/` directory **and the Prophet checkpoints**, so a project directory
needs only `input_automag.py`, the structure and that binary - and the machine
that runs it needs neither Python nor conda.

```bash
source /path/to/miniconda3/bin/activate automag
conda install -c conda-forge ninja -y
python -m pip install pyinstaller

./build.sh              # -> dist/run_automag
./build.sh --onedir     # -> dist/run_automag/ (starts faster)
./build.sh --clean      # remove build/ and dist/
```

`./build.sh` regenerates `build/deps.ninja` from the source tree (ninja has no
globbing) and runs `tools/build_executable.py`, which calls PyInstaller with
`--onefile`, the hidden imports of the dynamically loaded modules
(`pymatgen`, `ase`, `spglib`, `matplotlib`, `numpy`), `--collect-data` for
pymatgen/ase (their JSON/YAML tables are not covered by any PyInstaller hook)
and `--add-data` for the scripts of every stage directory, `common/`,
`calculators/`, `pipeline/` and `geometries/`.

Properties of the bundle worth knowing:

* only the `*.py` files of the stage directories are included - the example
  data some of them contain (`3_monte_carlo/Fe12O18/espins` alone is 600 MB of
  tracked ESpinS output) is skipped;
* `models/*.pt` (both Prophet checkpoints, 275 MB) are always added, and
  `common.paths.models_dir()` finds them at `<bundle>/_internal/models`;
* `prophet_spin` + torch are bundled **by default** (the frozen build then
  needs `--collect-all torch`, `--exclude-module torch._numpy` and
  `--collect-submodules ase`, otherwise torch's lazy imports and `ase.io.json`
  break).  `--no-prophet` drops them: the result is ~490 MB instead of ~5.7 GB;
* `external/codes.json` is copied in when it exists, so enumlib / VAMPIRE /
  ESpinS keep working from inside the bundle.

Without ninja, call PyInstaller directly:

```bash
python tools/build_executable.py               # one file (prophet included)
python tools/build_executable.py --onedir      # directory bundle, faster start
python tools/build_executable.py --no-prophet  # drop torch/prophet_spin (~490 MB)
python tools/build_executable.py --no-models   # prophet in, checkpoints outside
python tools/build_executable.py --dry-run     # print the command only
```

`--no-models` leaves the 275 MB of checkpoints out of the bundle; the binary then
reads them from `PROPHET_MODEL` / `PROPHET_HEAD` (or from an explicit
`model_path` in the stage params), which is convenient when a single `models/`
directory is shared over the network:

```bash
./build.sh --onedir                                    # ~5.7 GB, models inside
python tools/build_executable.py --onedir --name run_automag-extmodels --no-models
PROPHET_MODEL=/shared/models/prophet-v1-mag.pt \
PROPHET_HEAD=/shared/models/prophet-v1-mag-head.pt \
  ./dist/run_automag-extmodels/run_automag --dry-run
```

Use the binary like the Python entry point - no conda environment, no sources:

```bash
env -i PATH=/usr/bin:/bin ./dist/run_automag/run_automag --list-stages
env -i PATH=/usr/bin:/bin ./dist/run_automag/run_automag --dry-run --wait
AUTOMAG_EXE=/shared/bin/run_automag ./examples/vasp/run_pipeline_compiled.sh
```

A directory with only a structure file works too: the binary writes its default
`input_automag.py` and starts the Prophet magnetic search, using the models that
travel inside the bundle.

The example directories also ship `run_pipeline_slurm_compiled.sh`, which
submits the binary as a SLURM job; see `examples/README.md`.

Inside a frozen build `common.paths.code_root()` is the unpack directory, and
the runner starts stage scripts as `run_automag --run-script <path>` and
helpers as `--run-helper <path>` (relative to the bundle root, because a
one-file build unpacks into a fresh temporary directory on every start).
`external/codes.json` must be next to the binary, or be pointed at with
`AUTOMAG_EXTERNAL_CODES`.

Requirements of the build machine: the same Python you run automag with,
PyInstaller for that interpreter, and enough disk space in `build/` (a few
hundred MB).

---

### 5.1 What the bundle contains

| part | notes |
| --- | --- |
| Python + libraries | its own interpreter, no conda needed at run time |
| stage scripts | `0_conv_tests/`, `2_coll/`, ... without the example data |
| pymatgen / ase data | `periodic_table.json`, ... |
| Prophet checkpoints | `models/prophet-v1-mag.pt`, `models/prophet-v1-mag-head.pt` (always, unless `--no-models`) |
| `prophet_spin` + torch | bundled by default (`--no-prophet` leaves them out) |

Three ready-made variants:

| command | size | prophet |
| --- | --- | --- |
| `--onedir --name run_automag` | ~5.7 GB | inside, models included |
| `--onedir --name run_automag-extmodels --no-models` | ~5.4 GB | inside, models from `PROPHET_MODEL` / `PROPHET_HEAD` |
| `--onedir --name run_automag-noprophet --no-prophet --no-models` | ~500 MB | absent (DFT calculators only) |

Run the bundle on a machine without conda:

```bash
env -i PATH=/usr/bin:/bin ./dist/run_automag/run_automag --dry-run
```

## 6. External codes: enumlib, VAMPIRE, ESpinS

Three optional programs are used by the physics stages:

| code | used by | default command |
| --- | --- | --- |
| enumlib | `coll` (magnetic configurations) | `enum.x`, `makeStr.py` |
| VAMPIRE | `monte_carlo` (`mc_backend="vampire"`) | `vampire` |
| ESpinS | `monte_carlo` (`mc_backend="espins"`) | `mc.x` |

### Build them

```bash
sudo dnf install -y gcc gcc-c++ gfortran make git cmake   # or apt equivalents
./build.sh --external            # clones into external/src/, compiles, writes external/codes.json
./build.sh --all                 # the same + dist/run_automag
./build.sh --external --dry-run  # show the commands, do nothing
./build.sh --external --skip espins
```

Each code is cloned into `external/src/<name>` and built with its own
Makefile.  A code that fails to compile is recorded as `"available": false`
with a warning; the pipeline keeps working and only the step that needs that
code stops, telling the user what to do.  This is the expected outcome for
**ESpinS** with gfortran: its `make.sys` is written for Intel Fortran + MKL.

### Register codes you already have

Nothing is downloaded or compiled; the paths and the commands needed on your
cluster are written to `external/codes.json`:

```bash
./build.sh --external \
    --enumlib-path /opt/enumlib \
    --vampire-path /opt/vampire/vampire-serial \
    --vampire-pre-command  "module load intel" \
    --vampire-post-command "module purge" \
    --espins-path /opt/espins/mc.x
```

or keep the same information in a file and pass `--config cluster_codes.json`.

`external/codes.json` is read at run time by `common/external.py`;
`AUTOMAG_EXTERNAL_CODES` (file or directory) overrides its location, and
`PIPELINE["external_codes"]` overrides single entries per run.  Without any
`codes.json` the default commands above are used, i.e. the executables must be
on `PATH` - the historical behaviour.

---

## 7. Troubleshooting

| symptom | cause / fix |
| --- | --- |
| `TypeError: unsupported operand type(s) for \|` / `SyntaxError` on `X \| None` | Python is older than 3.10 - create a newer environment (§1, §2) |
| `ModuleNotFoundError: No module named 'pymatgen'` although it is installed | automag runs under a different interpreter than you installed into; stage steps inherit the interpreter of the process that runs `run_automag.py` - activate the environment first, or set `python` in the host section |
| `ModuleNotFoundError: No module named 'fireworks'` | the stage settings use `use_fireworks = True` (the historical default) - set `"submit": {"use_fireworks": False}` as all examples do |
| VASP stops immediately with a POTCAR error | `VASP_PP_PATH` is not set or lacks `potpaw_PBE` |
| `enum.x: command not found` | run `./build.sh --external`, register the path, or add it to `PATH` |
| VAMPIRE build takes forever | it is serial Fortran; expect ~40 min. Use `--vampire-path` with a pre-installed binary instead |
| ESpinS does not build | expected with gfortran: compile it with Intel Fortran + MKL and pass `--espins-path` |
| `error while loading shared libraries: libstdc++.so.6: version ...` | the code was compiled with the conda toolchain; use the system compilers (`tools/build_external.py` already prefers `/usr/bin`) |
| `ninja: not found` from `./build.sh` | `conda install -c conda-forge ninja`, `dnf install ninja-build`, or use `python tools/build_executable.py` |
| PyInstaller binary fails with a missing module | rebuild after installing that module (hidden imports are fixed at build time) |
| `ssh: Could not resolve hostname` / password prompt | password-less keys (`BatchMode=yes` is used) and the correct `host`/`user`/`identity_file` in `PIPELINE["hosts"]` |
| job is submitted but never appears | the `jobheader` must contain real `#SBATCH ...`/`#PBS ...` directives; automag only writes the script, the scheduler does the rest |

---

## 8. What to read next

* `README.md` - what the pipeline does, the input file, run modes, hosts,
  parallelism and every stage setting;
* `docs/AGENT_GUIDE.md` - the developer/agent reference (stage scripts,
  internal APIs, how to add a calculator);
* `examples/` - complete projects (`vasp`, `qe`, `fhi-aims`,
  `pymatgen_enumeration`, `prophet`, mixed-calculator setups and `remote_*`
  examples that run `2_coll` on another machine), each with the
  `run_pipeline*.sh` launchers: one command that waits (`run_pipeline.sh`),
  single passes for manual restarts (`run_pipeline_once.sh`,
  `run_pipeline_once_loop.sh`), the compiled variants and a self-submitting
  SLURM driver job.