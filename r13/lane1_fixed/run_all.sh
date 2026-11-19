#!/bin/zsh
# Re-run lane 1's I1-I8 groups from the patched copies (code-audit fixes + r-ball contact check). Sequential, nice 10.
cd research/tier2/union-closed/r13/lane1_fixed
echo "START $(date)" > logs/run_all.log
nice -n 10 <python> closure.py > logs/A_closure.log 2>&1; echo "A closure: $(grep -c PASS logs/A_closure.log) PASS, $(grep -c FAIL logs/A_closure.log) FAIL" >> logs/run_all.log
for g in SC K W Z C H G1 G3 G4 G2; do
  nice -n 10 <python> progA_checks.py $g > logs/A_$g.log 2>&1; echo "A $g: $(grep -c '^PASS' logs/A_$g.log) PASS, $(grep -c '^FAIL' logs/A_$g.log) FAIL $(date +%H:%M)" >> logs/run_all.log
done
for g in SC K W Z C H G1 G3 G4 G2; do
  nice -n 10 ./progB_fixed $g > logs/B_$g.log 2>&1; echo "B $g: $(grep -c '^PASS' logs/B_$g.log) PASS, $(grep -c '^FAIL' logs/B_$g.log) FAIL $(date +%H:%M)" >> logs/run_all.log
done
echo "END $(date)" >> logs/run_all.log
