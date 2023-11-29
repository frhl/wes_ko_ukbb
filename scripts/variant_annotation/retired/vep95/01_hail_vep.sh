#!/usr/bin/env bash
#
# @description annotate variants using Hail
# @depends quality controlled MatrixTables with variants.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=hail_vep
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/hail_vep.log
#SBATCH --error=logs/hail_vep.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-19
#SBATCH --requeue

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="data/unphased/wes/post-qc"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/vep/hail_hgvs/"

# hail script
readonly hail_script="scripts/01_hail_vep.py"

# input paths
readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} ) 
readonly in="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

if [ ! -f "${out_prefix}.ht" ]; then 
  set_up_hail
  set_up_vep
  set_up_pythonpath_legacy  
  python3 ${hail_script} \
       --chrom "${chr}" \
       --input_path "${in}" \
       --input_type "mt" \
       --out_prefix "${out_prefix}"
else
   raise_error "Hail VEP annotation for chr${chr} already exists. Skipping"
fi





