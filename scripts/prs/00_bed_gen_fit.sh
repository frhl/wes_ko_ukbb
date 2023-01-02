#!/usr/bin/env bash
#
# @description Generate plink files for indiviudals with whole exome sequencing data unavailable.
# These samples are used for training polygenic risks cores. 
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=bed_gen_fit
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/bed_gen_fit.log
#SBATCH --error=logs/bed_gen_fit.errors.log
#SBATCH --partition=long
#SBATCH --cpus-per-task 3
#SBATCH --array=1-22
#SBATCH --requeue

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/prs/00_bed_gen_fit.py"
readonly spark_dir="data/tmp/spark_dir"
readonly hap_dir="/well/lindgren/flassen/ressources/hapmap"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly out_dir="data/prs/hapmap/ukb_500k/fitting"
readonly out_prefix="${out_dir}/ukb_hapmap_500k_eur_maf5e-2_chr${chr}"
readonly hap_file="${hap_dir}/weights.l2.ldscore.liftover.ht"

readonly sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix}.bed" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  set -x
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --dataset "imp" \
     --ancestry "eur" \
     --hapmap ${hap_file} \
     --filter_missing 0.01 \
     --exclude_samples "${sample_list}" \
     --filter_to_unrelated_using_kinship_coef \
     --min_maf 0.01 \
     --liftover \
     --dbsnp \
     --out_prefix "${out_prefix}" \
     --out_type "plink" 
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




