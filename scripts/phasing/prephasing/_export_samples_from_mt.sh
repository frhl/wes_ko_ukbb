#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_samples_from_mt
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_samples_from_mt.log
#SBATCH --error=logs/export_samples_from_mt.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --array=21

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly hail_script="scripts/phasing/prephasing/_export_samples_from_mt.py"
readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/prephased/wes_union_calls/revision/50k"
readonly in_path="${in_dir}/ukb_shapeit5_whatshap_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/prephased/wes_union_calls/revision/50k"
readonly out_prefix="${out_dir}/ukb_shapeit5_whatshap_chr${chr}"

mkdir -p ${out_dir}

module purge
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --phased_path ${in_path} \
  --phased_type ${in_type} \
  --out_prefix ${out_prefix}

