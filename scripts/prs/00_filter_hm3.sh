#!/usr/bin/env bash
#
# @description Filter to samples with whole exome sequencing data available. These are the main
# samples that also have knockout status.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=filter_hm3
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/filter_hm3.log
#SBATCH --error=logs/filter_hm3.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-20

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/prs/00_filter_hm3.py"
readonly spark_dir="data/tmp/spark_dir"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly out_dir="data/imputed/hm3"
readonly out_prefix="${out_dir}/ukb_hm3_500k_chr${chr}"

readonly hap_dir="/well/lindgren/flassen/ressources/hapmap/ldpred2"
readonly hap_file="${hap_dir}/map_liftover.ht"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --hapmap ${hap_file} \
   --out_prefix "${out_prefix}" \
   --out_type "mt" \
   --chrom "${chr}"




