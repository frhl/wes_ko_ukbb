#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=sample_wes
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/sample_wes.log
#SBATCH --error=logs/sample_wes.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 3
#SBATCH --array=1,22

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/simulation/00_sample_wes.py"
readonly spark_dir="data/tmp/spark_dir"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="data/mt/annotated"
readonly in_file="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/simulation/mt"
readonly out_prefix="${out_dir}/ukb_unrel_samples"
readonly out_type="mt"

readonly seed="1995"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

# Sample WES
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --in_prefix "${in_file}"\
   --in_type "mt" \
   --filter_to_unrelated_using_kinship_coef \
   --out_prefix "${out_prefix}" \
   --out_type "${out_type}" 





