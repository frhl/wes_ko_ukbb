#!/usr/bin/env bash
#
# @description generate files of genotyped calls
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calls_gen
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calls_gen.log
#SBATCH --error=logs/calls_gen.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=20-22

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/02_geno_gen.py"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly out_dir="data/unphased/calls/new_liftover"
readonly out_prefix="${out_dir}/ukb_prefilter_calls_200k_chr${chr}"
readonly out_type="vcf"

readonly sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "$( get_hail_ext ${out_prefix} ${out_type})" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --chrom "${chr}" \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" \
     --extract_samples "${sample_list}" \
     --liftover \
     --min_mac 2 \
     --missing 0.05 \
     --ancestry "eur" \
     --dataset "calls"
else
  echo "file ${out} already exists. Skipping!"
fi


if [[ ${out_type} == "vcf" ]]; then
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix}.vcf.bgz" "tbi"
fi



