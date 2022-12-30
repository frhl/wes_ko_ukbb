#!/usr/bin/env bash
#
#$ -N submit_prs
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/submit_prs.log
#$ -e logs/submit_prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V



readonly script="scripts/prs/06_prs.sh"

qsub ${script}
sleep 15m 
qsub ${script}
sleep 15m
qsub ${script}
sleep 15m
qsub ${script}
sleep 15m
qsub ${script}
sleep 15m
qsub ${script}
sleep 15m
qsub ${script}






