#!/usr/bin/env bash
#
#$ -N test_hail
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/test_hail.log
#$ -e logs/test_hail.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q short.qe
#$ -t 22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="data/mt"
readonly vep_dir="data/vep/full/"
readonly spark_dir="data/tmp/spark"
readonly out_dir="derived/ptvs_test3"

# hail script
readonly hail_script="utils/hail_export.py"

# input path
readonly chr=${SGE_TASK_ID}
readonly in_phased="${in_dir}/ukb_wes_200k_annotated_chr${chr}.mt"
readonly in_unphased="${in_dir}/ukb_wes_200k_annotated_chr${chr}_singletons.mt"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_maf002_miss005_ptv_chr${chr}"
readonly out="${out_prefix}.mt"

# echo script into log
echo "##### SCRIPT ######" >> logs/test_hail.errors.log
echo $( cat ${hail_script} ) >> logs/test_hail.errors.log 
echo "##### END ######" >> logs/test_hail.errors.log

# run hail
set_up_hail
set_up_pythonpath
mkdir -p ${out_dir}
python3 "${hail_script}" 


print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





