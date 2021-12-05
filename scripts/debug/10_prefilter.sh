#!/usr/bin/env bash
#
#$ -N prefilter
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter.log
#$ -e logs/prefilter.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 24
#$ -q long.qc@@long.hga
#$ -t 21

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly in_dir="data/unphased/post-qc/"
readonly out_dir="data/phased/test-phasing"

readonly chr="${SGE_TASK_ID}"
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
readonly out_file="${out_dir}/ukb_wes_200k_refphased_chr${chr}.vcf.gz"




