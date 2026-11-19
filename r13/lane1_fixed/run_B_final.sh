#!/bin/zsh
cd ${0:A:h}
echo "START $(date)" > logs/run_B_final.log
for g in SC SCX K W Z ZX C H G1 G3 G4 G2; do
  nice -n 10 ./progB_final $g > logs/Bf_$g.log 2>&1; echo "Bf $g: $(grep -c '^PASS' logs/Bf_$g.log) PASS, $(grep -c '^FAIL' logs/Bf_$g.log) FAIL $(date +%H:%M)" >> logs/run_B_final.log
done
echo "END $(date)" >> logs/run_B_final.log
