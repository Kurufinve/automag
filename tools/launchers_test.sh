#!/bin/bash
# End-to-end check of the example launchers, without DFT.
#
#   A  run_pipeline.sh             one blocking command (--wait)
#   B  run_pipeline_once.sh        a single pass (--once)
#   C  run_pipeline_once_loop.sh   loop of passes, stops when the work is done
#
# The `conv_tests` stage is used because it is the cheapest one: two
# convergence points, no enumeration.  Only the submit phase is exercised --
# post-processing parses real VASP output, which a fake calculator cannot
# produce.  The calculator below copies a real OUTCAR so that the jobs are
# indistinguishable from finished VASP runs for the runner.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLE="${ROOT}/CalcFold/Al2Fe6Zr4/vasp/afm1/singlepoint/OUTCAR"
BASE=$(mktemp -d "${TMPDIR:-/tmp}/automag_launchers_XXXX")
FAIL=0

if [ ! -f "$SAMPLE" ]; then
    echo "SKIP: sample OUTCAR not found at $SAMPLE"
    exit 0
fi

FAKE="cp '${SAMPLE}' OUTCAR; cp POSCAR CONTCAR"

make_project() {  # target-dir
    cp -r "${ROOT}/examples/vasp/." "$1/"
    rm -rf "$1/__pycache__"
    python - "$1/input_automag.py" "$FAKE" <<'EOF'
import re, sys
path, fake = sys.argv[1], sys.argv[2]
text = open(path).read()
text = text.replace('"run": ["conv_tests", "coll", "relax", "monte_carlo", "mae"]',
                    '"run": ["conv_tests"]')
text = text.replace('"encut_values": range(500, 700, 50)', '"encut_values": [400]')
text = text.replace('"sigma_values": [0.05, 0.10]', '"sigma_values": [0.1]')
text = text.replace('"kpts_values": range(10, 40, 10)', '"kpts_values": [10]')
text = text.replace('"calculator_command": "module load vasp/6.4.3; mpirun vasp_std"',
                    '"calculator_command": %r' % fake)
text = re.sub(r'"jobheader":\s*"""..*?"""', '"jobheader": "#!/bin/bash"', text, flags=re.S)
open(path, "w").write(text)
EOF
}

check() {  # name exit-code dir
    local name=$1 status=$2 dir=$3 incar outcar
    incar=$(find "$dir/CalcFold" -name INCAR 2>/dev/null | wc -l)
    outcar=$(find "$dir/CalcFold" -name OUTCAR 2>/dev/null | wc -l)
    if [ "$status" -eq 0 ] && [ "$incar" -ge 1 ] && [ "$outcar" -ge 1 ]; then
        echo "PASS ${name}: exit 0, ${incar} input(s), ${outcar} result(s)"
    else
        echo "FAIL ${name}: exit ${status}, ${incar} input(s), ${outcar} result(s)"
        tail -8 "$dir/run.log"
        FAIL=1
    fi
}

run_case() {  # name launcher extra-env
    local name=$1 launcher=$2 dir="${BASE}/$3"
    make_project "$dir"
    ( cd "$dir" && AUTOMAG_CONDA_ENV="" AUTOMAG_ONCE_INTERVAL=1 \
        AUTOMAG_ONCE_MAX_ROUNDS=4 AUTOMAG_HOME="$ROOT" \
        timeout 900 "./${launcher}" --stages conv_tests --steps submit \
        > run.log 2>&1 )
    check "$name" $? "$dir"
}

run_case "run_pipeline.sh (--wait)"          run_pipeline.sh           blocking
run_case "run_pipeline_once.sh (--once)"     run_pipeline_once.sh      once
run_case "run_pipeline_once_loop.sh (loop)"  run_pipeline_once_loop.sh loop

echo "projects: ${BASE}"
exit "$FAIL"