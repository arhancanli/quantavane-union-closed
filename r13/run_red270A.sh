#!/bin/zsh
cd research/tier2/union-closed/r13
for f in eps280:N5 noPhi:N4 noK2:N3 theta32:N3 Psi_low:N3; do n=${f%%:*}; w=${f##*:}; echo "== red270 $n ($w)"; nice -n 10 <python> certA_jk.py red270_$n.json $w 2>&1 | tail -1; done
echo DONE
