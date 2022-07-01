#!/usr/bin/env bash
#
#$ -N phase_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/phase_chunks.log
#$ -e logs/phase_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 19
#$ -V

set -o errexit
set -o nounset

#source utils/qsub_utils.sh
source utils/hail_utils.sh
#source utils/vcf_utils.sh

set +eu
set_up_hail
#set_up_pythonpath_legacy
set -eu


