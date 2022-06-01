#!/usr/bin/env bash
#
#
#$ -N merge_common_rare_markers
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/merge_common_rare_markers.log
#$ -e logs/merge_common_rare_markers.errors.log
#$ -P lindgren.prjc
#$ -q test.qc
#$ -pe shmem 1
#$ -t 1-22
#$ -V


set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/combined/01_combine_ko_rare_common.py"

readonly chr="${SGE_TASK_ID}"
readonly variants_dir="data/mt/annotated"

# note: assuming ko and rare variants have already been merged
readonly ko_rare_dir="data/conditional/rare/combined"
readonly common_dir="data/conditional/common/marker_mt"
readonly out_dir="data/conditional/combined"

readonly ko_rare_path="${ko_rare_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense_markers.txt.gz"
readonly common_path="${common_dir}/conditional_markers_chrundefined_markers.txt.gz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_markers.txt"

# make header
zcat $ko_rare_path | head -n1 > ${out_prefix}

# append files and subset to right chromosome
zcat $common_path $ko_rare_path | grep -v locus | grep chr${chr} >> ${out_prefix}



