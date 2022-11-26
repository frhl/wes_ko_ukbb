#!/usr/bin/env bash
#
#$ -N tally_dir
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o tally.txt
#$ -e tally.txt
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

readonly thetime=`date`
echo "### ${thetime} ###"
echo $( du -sh data/* | sort -hr )
echo "### (end) ###"




