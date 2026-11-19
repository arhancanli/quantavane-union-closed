#!/bin/zsh
cd research/tier2/union-closed/r13
nice -n 10 <python> certA_jk.py params_e270.json N5 N1 N2 N4 N3 > certA_e270.log 2>&1
nice -n 10 <python> certB_jk.py params_e270.json N5 N1 N2 N4 N3 > certB_e270.log 2>&1
echo DONE >> certB_e270.log
