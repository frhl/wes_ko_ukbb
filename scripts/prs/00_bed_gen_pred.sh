#!/usr/bin/env bash
#
# @description Filter to samples with whole exome sequencing data available. These are the main
# samples that also have knockout status.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=bed_gen_pred
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/bed_gen_pred.log
#SBATCH --error=logs/bed_gen_pred.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/prs/00_bed_gen_pred.py"
readonly spark_dir="data/tmp/spark_dir"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly out_dir="data/prs/hapmap/ukb_500k_hm3/validation"
readonly out_prefix="${out_dir}/ukb_hapmap_500k_eur_chr${chr}"

readonly common_dir="data/unphased/imputed/liftover"
readonly common_path="${common_dir}/ukb_imp_200k_chr${chr}.mt"

readonly hap_dir="/well/lindgren/flassen/ressources/hapmap/ldpred2"
readonly hap_file="${hap_dir}/map_liftover.ht"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --hapmap ${hap_file} \
     --common_path "${common_path}" \
     --min_maf 0.01 \
     --filter_missing 0.05 \
     --filter_invariant \
     --out_prefix "${out_prefix}" \
     --out_type "plink" 
else
  print_update "file ${out} already exists. Skipping!"
fi




