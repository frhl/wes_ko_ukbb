#!/usr/bin/env bash
#
# script that counts number of heterozygote alleles in a gene
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=count_hets
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/count_hets.log
#SBATCH --error=logs/count_hets.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1-22
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/permute/02_count_hets.py"

readonly in_dir="data/mt/annotated"
readonly out_dir="data/permute/counts"

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly input_path="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.mt"
readonly input_type='mt'

readonly csqs="pLoF,damaging_missense"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_phased_counts_chr${chr}"
readonly out_type="mt"

mkdir -p ${out_dir}

if [ ! -d "${out_prefix}.mt" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 ${hail_script} \
    --input_path ${input_path} \
    --input_type ${input_type} \
    --out_prefix ${out_prefix} \
    --out_type ${out_type} \
    --csqs_category ${csqs} \
    --discard_unphased \
    --use_loftee
else
  >&2 echo "${out_prefix}.mt already exists. Skipping.."
fi




