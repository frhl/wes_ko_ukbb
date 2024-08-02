#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=liftover_hm3
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/liftover_hm3.log
#SBATCH --error=logs/liftover_hm3.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/prs/00_liftover_hm3.py"
readonly spark_dir="data/tmp/spark_dir"

readonly hap_dir="/well/lindgren/flassen/ressources/hapmap/ldpred2"
readonly hm3_file="${hap_dir}/map.txt.gz"
readonly hm3_out_prefix="${hap_dir}/map_liftover"
readonly hm3plus_file="${hap_dir}/map_hm3_plus.txt.gz"
readonly hm3plus_out_prefix="${hap_dir}/map_hm3_plus_liftover"

mkdir -p ${spark_dir}
mkdir -p ${hap_dir}

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --hapmap ${hm3_file} \
   --out_prefix "${hm3_out_prefix}"

python3 "${hail_script}" \
   --hapmap ${hm3plus_file} \
   --out_prefix "${hm3plus_out_prefix}"





