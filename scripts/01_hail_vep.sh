#!/usr/bin/env bash
#
#$ -N hail_vep
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/hail_vep.log
#$ -e logs/hail_vep.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qc
#$ -t 23

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="data/unphased/wes/post-qc"
readonly spark_dir="data/tmp/spark"
readonly vep_dir="data/vep/full"
readonly out_dir="data/vep/hail"

# hail script
readonly hail_script="scripts/01_hail_vep.py"

# input paths
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

if [ ! -f "${out_prefix}.ht" ]; then 
  SECONDS=0
  set_up_hail
  set_up_vep
  set_up_pythonpath_legacy  
  python3 "${hail_script}" \
       --chrom "${chr}" \
       --input_path "${in}" \
       --input_type "mt" \
       --vep_path "${vep}" \
       --out_prefix "${out_prefix}" \
       && print_update "Finished Hail VEP annotation chr${chr}" ${SECONDS} \
       || raise_error "Hail VEP annotation for chr${chr} failed"
else
   raise_error "Hail VEP annotation for chr${chr} already exists. Skipping"
fi





