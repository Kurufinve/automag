#!/bin/bash
# Build automag.
#
#   ./build.sh                  # single-file executable -> dist/run_automag
#   ./build.sh --external       # enumlib / ESpinS / VAMPIRE -> external/
#   ./build.sh --all            # both
#   ./build.sh --onedir         # directory bundle (faster start)
#   ./build.sh --clean          # remove build/, dist/
#
# External codes (see tools/build_external.py):
#
#   ./build.sh --external                       # download + compile
#   ./build.sh --external --vampire-path /opt/vampire/vampire-serial \
#                          --espins-path  /opt/espins/mc.x
#        # nothing is downloaded; the paths (and optional
#        # --<code>-pre-command / --<code>-post-command) are written to
#        # external/codes.json, which automag reads at run time.
#
# A code that fails to compile is recorded as unavailable with a warning:
# compile it manually and rerun with --<code>-path.
#
# Requires ninja and pyinstaller in the active conda environment:
#   source /path/to/miniconda3/bin/activate automag
#   conda install -c conda-forge ninja && pip install pyinstaller
set -euo pipefail

cd "$(dirname "$0")"

case "${1:-}" in
    --clean)
        rm -rf build dist
        echo "Removed build/ and dist/ (external/ is kept; remove it by hand)."
        exit 0
        ;;
    --external)
        shift
        exec python tools/build_external.py "$@"
        ;;
    --all)
        shift
        python tools/build_external.py "$@"
        exec python tools/build_executable.py
        ;;
    --onedir)
        exec python tools/build_executable.py --onedir
        ;;
esac

if ! command -v ninja >/dev/null 2>&1; then
    echo "ERROR: ninja not found. Install it with:" >&2
    echo "    conda install -c conda-forge ninja" >&2
    echo "(or run './build.sh --external' / 'python tools/build_executable.py')" >&2
    exit 1
fi

# Regenerate the source list before the first build so that `include
# build/deps.ninja` never refers to a missing file.
mkdir -p build
python tools/write_ninja_deps.py --output build/deps.ninja

exec ninja "$@"
