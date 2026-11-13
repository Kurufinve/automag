#!/bin/bash
# End-to-end check of the "remote" examples without a second machine.
#
# The examples describe a cluster reachable by `ssh`; here the same host
# definition is switched to the `shell` transport (a wrapper that runs the
# command locally), which exercises the identical code path in
# pipeline/hosts.py: everything is wrapped, cd'ed into the host workdir and
# executed through the wrapper.  A fake calculator writes outputs instead of
# running DFT.
#
# Usage: bash tools/remote_examples_test.sh [example ...]
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BASE=$(mktemp -d "${TMPDIR:-/tmp}/automag_remote_XXXX")
FAIL=0

FAKE_VASP="sh -c 'printf \"General timing and accounting informations for this job:\\\\n\" > OUTCAR; printf \"       N       E                     dE             d eps       ncg     rms          rms(c)\\\\nDAV:   1   -400.0   -400.0   -0.456E+03   1234   0.123E+00\\\\n   1 F= -400.0 E0= -400.0  d E =0.0  mag= 4.0 -4.0 4.0 -4.0 0.0 0.0 0.0 0.0 0.0 0.0\\\\n\" > OSZICAR'"

for example in "${@:-remote_vasp remote_fhi_aims remote_qe}"; do
    for name in $example; do
        src="${ROOT}/examples/${name}"
        if [ ! -d "$src" ]; then
            echo "FAIL ${name}: no such example"; FAIL=1; continue
        fi
        dir="${BASE}/${name}"
        cp -r "$src/." "$dir/"
        rm -rf "$dir/__pycache__"
        python - "$dir/input_automag.py" "$FAKE_VASP" <<'EOF'
import re, sys
path, fake = sys.argv[1], sys.argv[2]
text = open(path).read()
# ssh -> shell transport through a wrapper, scheduler off (no real cluster)
text = re.sub(
    r'("cluster": \{\n\s*)"transport": "ssh",.*?\n(\s*\},)',
    r'\1"transport": "shell",\n            "wrapper": "sh -c",\n'
    r'            "workdir": "%s",\n            "scheduler": "none",\n\2' % __import__("os").path.dirname(path),
    text, count=1, flags=re.S)
# no batch directives: the calculations run as plain processes on the host
text = re.sub(r'"jobheader":\s*"""..*?"""', '"jobheader": "#!/bin/bash"', text,
              flags=re.S)
# a calculator that writes its outputs instead of running DFT
if 'examples' in path or True:
    text = re.sub(r'"calculator_command":\s*"[^"]*"',
                  '"calculator_command": %r' % fake, text, count=1)
open(path, "w").write(text)
EOF
        ( cd "$dir" && AUTOMAG_HOME="$ROOT" AUTOMAG_CONDA_ENV="" \
            timeout 900 python "$ROOT/run_automag.py" --execution process --wait \
            > run.log 2>&1 )
        status=$?
        # the plan/running lines name the host of the step
        host=$(grep -oE "@[a-z_]+" "$dir/run.log" | sort -u | tr -d '@' | paste -sd, -)
        out=$(find "$dir/CalcFold" -name "OUTCAR" 2>/dev/null | wc -l)
        if [ "$status" -eq 0 ] && [ "$host" = "cluster" ]; then
            echo "PASS ${name}: exit 0, ran on host '${host}', ${out} output(s)"
        else
            echo "FAIL ${name}: exit ${status}, hosts '${host}'"
            tail -12 "$dir/run.log"
            FAIL=1
        fi
    done
done

echo "projects: ${BASE}"
exit "$FAIL"
