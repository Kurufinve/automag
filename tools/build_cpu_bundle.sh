#!/usr/bin/env bash
# Build the small (CPU-only) automag bundles.
#
# The 5.7 GB bundle is big because the CUDA build of torch brings ~3 GB of GPU
# libraries that CPU-only machines never use.  A CPU wheel of the *same* torch
# version gives identical results and fits in ~350 MB, so the bundle becomes
# ~830 MB once the checkpoints are kept outside and torch's own test binaries
# are dropped.
#
# Usage:
#   tools/build_cpu_bundle.sh                # both bundles, models inside + outside
#   tools/build_cpu_bundle.sh --no-models    # only the small one (models outside)
#   tools/build_cpu_bundle.sh --setup-only   # just prepare the CPU environment
#
# Environment:
#   AUTOMAG_CPU_ENV      venv to build from (default ~/.cache/automag/cpuenv)
#   AUTOMAG_CPU_PYTHON   interpreter used to create that venv (default python3)
#   AUTOMAG_TORCH_CPU    torch version to install (default 2.14.1)
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
ENV_DIR=${AUTOMAG_CPU_ENV:-$HOME/.cache/automag/cpuenv}
# The bundle must match the interpreter of the reference environment: prophet
# wheels are built for one ABI, so prefer 3.12 over whatever python3 is.
if [ -n "${AUTOMAG_CPU_PYTHON:-}" ]; then
    BASE_PYTHON=$AUTOMAG_CPU_PYTHON
elif command -v python3.12 >/dev/null 2>&1; then
    BASE_PYTHON=$(command -v python3.12)
else
    BASE_PYTHON=$(command -v python3)
fi
TORCH_VERSION=${AUTOMAG_TORCH_CPU:-2.14.1}

BUILD_EMBEDDED=1
SETUP_ONLY=0
for argument in "$@"; do
    case "$argument" in
        --no-models) BUILD_EMBEDDED=0 ;;
        --setup-only) SETUP_ONLY=1 ;;
        *) echo "unknown option: $argument" >&2; exit 2 ;;
    esac
done

if [ ! -x "$ENV_DIR/bin/python" ]; then
    echo "=== creating the CPU build environment: $ENV_DIR"
    "$BASE_PYTHON" -m venv "$ENV_DIR"
fi
# The tag comes from the environment that will actually run the build.
tag=$("$ENV_DIR/bin/python" -c 'import sys; print("cp%d%d" % sys.version_info[:2])')
echo "=== build environment: $ENV_DIR ($tag)"
PIP="$ENV_DIR/bin/pip"
"$PIP" install -q --upgrade pip

echo "=== installing torch $TORCH_VERSION+cpu"
if ! "$PIP" install -q --force-reinstall "torch==$TORCH_VERSION+cpu" \
        --index-url https://download.pytorch.org/whl/cpu; then
    url_name="torch-${TORCH_VERSION}%2Bcpu-${tag}-${tag}-manylinux_2_28_x86_64.whl"
    file_name="torch-${TORCH_VERSION}+cpu-${tag}-${tag}-manylinux_2_28_x86_64.whl"
    cache=$ENV_DIR/wheels
    mkdir -p "$cache"
    if [ ! -f "$cache/$file_name" ]; then
        echo "    the CPU index has no $TORCH_VERSION+cpu for $tag, downloading the wheel"
        curl -sSL --max-time 900 -o "$cache/$file_name" \
            "https://download-r2.pytorch.org/whl/cpu/$url_name"
    fi
    "$PIP" install -q --force-reinstall --no-deps "$cache/$file_name"
fi

echo "=== installing the project requirements and the build tools"
"$PIP" install -q -r "$ROOT/requirements.txt" pyinstaller ninja e3nn matscipy torch_geometric

# prophet-mlip is distributed as a wheel inside the Prophet repository, not from
# PyPI: reuse the copy that the reference environment already has.
source_env=${AUTOMAG_SOURCE_ENV:-$(python -c 'import sys, os; print(os.path.dirname(os.path.dirname(sys.executable)))' 2>/dev/null || true)}
if [ -n "${source_env:-}" ] && [ -d "$source_env/lib" ]; then
    site=$("$ENV_DIR/bin/python" -c 'import site; print(site.getsitepackages()[0])')
    for package in prophet prophet_spin; do
        if [ -d "$source_env/lib/python3."*/site-packages/$package ] && \
           [ ! -d "$site/$package" ]; then
            cp -a "$source_env"/lib/python3.*/site-packages/$package "$site/"
            echo "    copied $package from $source_env"
        fi
    done
    for info in "$source_env"/lib/python3.*/site-packages/prophet_mlip-*.dist-info; do
        [ -d "$info" ] && cp -a "$info" "$site/" 2>/dev/null || true
    done
fi

"$ENV_DIR/bin/python" - <<'PY'
import torch
print(f"    torch {torch.__version__}, cuda available: {torch.cuda.is_available()}")
import prophet_spin  # noqa: F401 - the build needs it importable
print("    prophet_spin importable")
PY

if [ "$SETUP_ONLY" = "1" ]; then
    echo "=== environment ready: $ENV_DIR"
    exit 0
fi

cd "$ROOT"
if [ "$BUILD_EMBEDDED" = "1" ]; then
    echo "=== building dist/run_automag-cpu (checkpoints inside)"
    rm -rf dist/run_automag-cpu
    "$ENV_DIR/bin/python" tools/build_executable.py --onedir --name run_automag-cpu
fi

echo "=== building dist/run_automag-cpu-extmodels (checkpoints outside)"
rm -rf dist/run_automag-cpu-extmodels
"$ENV_DIR/bin/python" tools/build_executable.py --onedir \
    --name run_automag-cpu-extmodels --no-models

du -sh dist/run_automag-cpu* 2>/dev/null || true
cat <<'EOF'

Use the small bundle with:

    export PROPHET_MODEL=/path/to/models/prophet-v1-mag.pt
    export PROPHET_HEAD=/path/to/models/prophet-v1-mag-head.pt
    AUTOMAG_EXE=<repo>/dist/run_automag-cpu-extmodels/run_automag-cpu-extmodels \
        ./run_pipeline_compiled.sh

Or, without bundling torch at all, keep the 359 MB bundle and point
AUTOMAG_HELPER_PYTHON at any environment that has torch + prophet_spin.
EOF
