#!/usr/bin/env bash
#
# @description merge phased vcf with unphased parents for later switch error calculation.
# @note - takes about ~ 24h with 1 a core
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=extract_prephase_sites
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/extract_prephase_sites.log
#SBATCH --error=logs/extract_prephase_sites.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 5
#SBATCH --array=20,22

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly hail_script="scripts/phasing/03_extract_sites.py"
readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly prephased_dir="data/prephased/wes_union_calls"
readonly prephased_path="${prephased_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.gz"
readonly prephased_type="vcf"

readonly out_dir="data/prephased/wes_union_calls"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"

mkdir -p ${out_dir}


# note: we can't combine data using bcftools beacuse some variants
# share the same position (but not same reference alleles)
module purge
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --prephased_path ${prephased_path} \
  --prephased_type ${prephased_type} \
  --out_prefix ${out_prefix}

