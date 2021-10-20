#!/usr/bin/env bash
#
#$ -N debug_hail_vep
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/debug_hail_vep.log
#$ -e logs/debug_hail_vep.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q short.qf
#$ -t 21

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="/gpfs3/well/lindgren/UKBIOBANK/nbaya/wes_200k/phase_ukb_wes/data/old-all/all_variants/unphased"
readonly spark_dir="data/tmp/spark"
readonly vep_dir="data/vep/full"
readonly out_dir="data/debug_tmp"

# hail script
readonly hail_script="scripts/01_hail_vep.py"

# input paths
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in="${in_dir}/old_ukb_wes_hail_chr${chr}.vcf.bgz"
readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

# run hail
set_up_hail
set_up_vep
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --chrom ${chr} \
     --input_path ${in}\
     --input_type "vcf" \
     --vep_path ${vep} \
     --out_prefix ${out_prefix}

print_update "Finished running HAIL for chr${chr}" "${SECONDS}"





