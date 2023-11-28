#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=vep105
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/vep.log
#SBATCH --error=logs/vep.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-21

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/unphased/wes/post-qc"
readonly in="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"

readonly json_path="utils/configs/vep105_revel_float.json"

readonly out_dir="data/vep/vep105/vep_out"
readonly out_prefix="${out_dir}/UKB.chr${chr}.exome_array.variants_only.vep105"
readonly hail_script="scripts/variant_annotation/vep105/01_vep.py"

mkdir -p ${out_dir}
mkdir -p ${spark_dir}

if [ ! -f "${out_prefix}_vep.ht/_SUCCESS" ]; then
  set_up_hail 0.2.97
  set_up_vep105
  set_up_pythonpath_legacy
  python3 ${hail_script} \
       --input_path "${in}" \
       --out_prefix "${out_prefix}" \
       --json_path "${json_path}"
else
  >&2 echo "${out_prefix}* already exists."
fi

