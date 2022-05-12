#!/usr/bin/env bash
#
#
#$ -N _test_conda
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_test_conda.log
#$ -e logs/_test_conda.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -t 1-40
#$ -q test.qc

set -o errexit
set -o nounset

source utils/bash_utils.sh

set_up_RSAIGE 1.0.4


