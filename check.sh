#!/bin/sh
# One command:
#   sh check.sh          standard check: Program A on every finite inequality N1-N6 (about 15-30 minutes), Program B on
#                        N5, N1, N2, N6, the red control for N5, Proposition P, and the translate computations (small groups)
#   sh check.sh full     also Program B on N4 and N3 (about 3 hours)
# Needs one Python with python-flint (0.9.0 used here) and mpmath (1.3.0 used here): set PYTHON=/path/to/python.
# The certified inputs (I1)-(I9) and Boppana's inequality (I6) have their own programs and logs (see README.md).
set -eu
cd "$(dirname "$0")"
MODE=${1:-standard}
PY=${PYTHON:-python3}
fail() { echo "CHECK FAILED: $*"; exit 1; }
"$PY" -c 'import flint, mpmath' 2>/dev/null || fail "PYTHON needs python-flint and mpmath (pip install python-flint==0.9.0 mpmath)"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

echo "1. Program A (python-flint, Arb balls), parameters r13/params_e290.json (eps = 29/100000): N5 N1 N2 N6 N4 N3"
"$PY" -I r13/certA_cut.py r13/params_e290.json N5 N1 N2 N6 N4 N3 | tee "$T/A.txt" | grep -E "^N[0-9]:"
for g in N5 N1 N2 N6 N4 N3; do grep -q "^$g: PASS" "$T/A.txt" || fail "Program A, $g"; done

echo "2. Program B (binary64 outward intervals + mpmath.iv), written separately"
if [ "$MODE" = full ]; then GB="N5 N1 N2 N6 N4 N3"; else GB="N5 N1 N2 N6"; fi
"$PY" -I r13/certB_cut.py r13/params_e290.json $GB | tee "$T/B.txt" | grep -E "^N[0-9]:"
for g in $GB; do grep -q "^$g: PASS" "$T/B.txt" || fail "Program B, $g"; done

echo "3. Red control: eps = 59/200000 (2.95e-4) must fail N5 in both programs"
"$PY" -I r13/certA_cut.py r13/red285_eps295.json N5 | grep -q "^N5: FAIL" || fail "Program A accepted eps = 2.95e-4"
"$PY" -I r13/certB_cut.py r13/red285_eps295.json N5 | grep -q "^N5: FAIL" || fail "Program B accepted eps = 2.95e-4"
echo "   refused by both"

echo "4. Proposition P (barrier): PASS at eps = 1205/10^6, FAIL just below eps_P (12044/10^7), in both programs"
"$PY" -I barrier/d2_certP_A.py 1205 1000000 | grep -q "RESULT: PASS" || fail "Proposition P, program A"
"$PY" -I barrier/d2_certP_B.py 1205 1000000 | grep -q "RESULT: PASS" || fail "Proposition P, program B"
"$PY" -I barrier/d2_certP_A.py 12044 10000000 | grep -q "RESULT: FAIL" || fail "Proposition P red control, program A"
"$PY" -I barrier/d2_certP_B.py 12044 10000000 | grep -q "RESULT: FAIL" || fail "Proposition P red control, program B"
"$PY" -I barrier/referee/p_barrier.py > "$T/P.txt" 2>&1 || fail "Proposition P, referee program"
grep -q "eps in \[1205e-6, m-1/2\]): PASS" "$T/P.txt" || fail "Proposition P, referee range check"
echo "   ok"

echo "5. Translates: exhaustive check of the reciprocity and of Frankl's condition in Z6, S3, D4, Q8, D5"
"$PY" -I translates/referee/groups_check.py small > "$T/G.txt"
[ "$(grep -c 'step failures: 0 \[\]' "$T/G.txt")" = 5 ] || fail "translate step failures"
[ "$(grep -c 'Frankl failures: 0' "$T/G.txt")" = 5 ] || fail "translate Frankl failures"
echo "   ok"
echo "ALL CHECKS PASSED ($MODE): c_UC >= c* + 29/100000 = 0.38317525996679779...; Proposition P; translate computations."
