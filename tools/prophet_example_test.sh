#!/bin/bash
# End-to-end check of examples/prophet: the whole coll stage (enumeration ->
# prophet_spin energies -> postprocess) with the real Prophet-V1-Mag model on
# CPU, then the stage results file.
set -uo pipefail

ROOT=/home/dpoletaev/magnetic/automag
PROJ=$(mktemp -d "${TMPDIR:-/tmp}/automag_prophet_XXXX")
cp -r "$ROOT/examples/prophet/." "$PROJ/"
rm -rf "$PROJ/__pycache__"

python - "$PROJ/input_automag.py" <<'EOF'
import re, sys
path = sys.argv[1]
text = open(path).read()
# no batch directives: run the calculations as plain local processes
text = re.sub(r'"jobheader":\s*"""..*?"""', '"jobheader": "#!/bin/bash"', text, flags=re.S)
# the example's environment_activate is a placeholder
text = re.sub(r'"environment_activate":\s*\(.*?\),',
              '"environment_activate": "",', text, flags=re.S)
text = text.replace('"environment_deactivate": "conda deactivate"',
                    '"environment_deactivate": ""')
text = text.replace('"run": ["coll", "monte_carlo"]', '"run": ["coll"]')
open(path, "w").write(text)
EOF

cd "$PROJ" || exit 1
echo "project: $PROJ"
echo "model:   ${PROPHET_MODEL:-/opt/uiap/prophet/prophet-v1-mag.pt}"
AUTOMAG_HOME="$ROOT" AUTOMAG_CONDA_ENV="" PROPHET_DEVICE="${PROPHET_DEVICE:-cpu}" \
    timeout 3000 python "$ROOT/run_automag.py" --wait > run.log 2>&1
status=$?
echo "exit=$status"
tail -12 run.log

table=$(ls "$PROJ"/CalcFold/*singlepoint*.txt 2>/dev/null | head -1)
echo "results table: ${table:-MISSING}"
if [ -n "$table" ]; then
    cat "$table"
fi
done_count=$(grep -c "done" "$table" 2>/dev/null || true)
echo "---"
if [ "$status" -eq 0 ] && [ -n "$table" ] && grep -q "energy=" "$table"; then
    echo "PROPHET E2E OK"
    exit 0
fi
echo "PROPHET E2E FAILED"
tail -25 "$PROJ"/2_coll/.automag_logs/*.log 2>/dev/null
exit 1