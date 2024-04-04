#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=extract_marker_gt
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/extract_marker_gt.log
#SBATCH --error=logs/extract_marker_gt.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 5

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/common/08_extract_marker_gt.py"

readonly out_dir="data/conditional/common/markers/2024"
readonly out_prefix="${out_dir}/common_conditional"
readonly out_checkpoint="${out_prefix}_checkpoint.mt"
readonly out_type="vcf"

readonly markers_dir="data/conditional/common/spa_iter_new_2024"
readonly markers=$(cat ${markers_dir}/*.markers | cut -f2 | sort -u | tr "\n" ",")
readonly chroms=$(cat ${markers_dir}/*.markers | sed '/^[[:space:]]*$/d' | cut -d":" -f1)

mkdir -p ${out_dir}

# create new file containing all the aggregated results
cat ${markers_dir}/*.markers > "${out_prefix}.markers"


# Get matrix filtered to common variants that pass our conditional analysis
if [ ! -d "${out_prefix}.mt" ]; then 
  SECONDS=0
  set +u
  set_up_hail
  set -u
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --markers ${markers} \
     --out_type ${out_type} \
     --out_prefix ${out_prefix} \
     && print_update "Finished filtering imputed genotypes ${out_prefix}" ${SECONDS} \
     || raise_error "Filtering imputed genotypes for for ${out_prefix} failed!"
else
  >&2 echo "${out_prefix} already exists. Skipping"
fi




