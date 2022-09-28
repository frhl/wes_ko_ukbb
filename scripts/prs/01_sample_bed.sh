#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=sample_bed
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/sample_bed.log
#SBATCH --error=logs/sampled_bed.errors.log
#SBATCH --partition=long
#SBATCH --cpus-per-task 3
#SBATCH --array=1-22
#SBATCH --requeue

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/prs/00_bed_gen.py"
readonly spark_dir="data/tmp/spark_dir"
readonly hap_dir="/well/lindgren/flassen/ressources/hapmap"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly out_dir="data/prs/hapmap/ld/unrel_kin_eur_10k"
readonly out_prefix="${out_dir}/short_ukb_hapmap_rand_10k_eur_chr${chr}"
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
     --exclude_related \
     --filter_missing 0.01 \
     --random_samples 10000 \
     --filter_to_unrelated_using_kinship_coef \
     --exclude_samples "${sample_list}" \
     --write_samples \
     --hapmap "${hap_file}" \
     --min_maf 0.01 \
     --liftover \
     --dbsnp \
     --only_valid_contigs \
     --out_prefix "${out_prefix}" \
     --out_type "plink" 
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




