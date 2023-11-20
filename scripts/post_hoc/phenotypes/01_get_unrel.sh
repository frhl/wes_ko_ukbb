#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_unrel
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_unrel.log
#SBATCH --error=logs/get_unrel.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/post_hoc/01_get_unrel.py"
readonly spark_dir="data/tmp/spark_dir"

readonly in_dir="data/knockouts/alt/pp90/combined"
readonly in_file="${in_dir}/ukb_eur_wes_200k_chr21_pLoF_damaging_missense.vcf.bgz"
readonly in_type="vcf"

readonly out_dir="data/post_hoc/unrelated"
readonly out_prefix="${out_dir}/ukb_wes_ko_samples"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --in_file ${in_file} \
   --in_type ${in_type} \
   --out_prefix "${out_prefix}"




