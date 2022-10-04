#!/usr/bin/env bash
#
# @description Prune phased chunk boundaries so that only a few overlapping sites remain.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=trim_chunks
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/trim_chunks.log
#SBATCH --error=logs/trim_chunks.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 6
#SBATCH --array=2

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/phasing/05_trim_chunks.py"
readonly spark_dir="data/tmp/spark"


readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly main_dir="data/phased/wes_union_calls/chunks/shapeit4"
readonly in_dir="${main_dir}/ukb_eur_wes_union_calls_200k_chr${chr}-16xshort/"
readonly in_prefix="shapeit4_prs100000_pro50000_mprs100000"

readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"

readonly out_dir="data/phased/wes_union_calls/trimmed"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly out="${out_prefix}.vcf.bgz"
readonly trio="${out_prefix}.trio"

mkdir -p ${out_dir}

if [ ! -f ${out} ]; then
  set_up_hail
  set_up_pythonpath_legacy
  SECONDS=0
  set -x
  python3 ${hail_script} \
      --in_dir "${in_dir}" \
      --in_ext ".vcf.gz" \
      --in_prefix "${in_prefix}" \
      --new_overlap_size 5000 \
      --out_prefix "${out_prefix}" \
      --out_type "vcf" \
      && print_update "Finished merging phased data for chr${chr}" ${SECONDS} \
      || raise_error "Merging phased data for chr${chr} failed" 
  set +x
else
    print_update "Warning: ${out} already exists! Skipping." | tee /dev/stderr
fi






