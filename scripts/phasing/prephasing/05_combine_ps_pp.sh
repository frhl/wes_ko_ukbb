#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=combine_ps_pp
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/combine_ps_pp.log
#SBATCH --error=logs/combine_ps_pp.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 3
#SBATCH --array=20,22
#
#$ -N combine_ps_pp
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/combine_ps_pp.log
#$ -e logs/combine_ps_pp.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 20-22
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly hail_script="scripts/phasing/prephasing/05_combine_ps_pp.py"
readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly ref_dir="data/prephased/wes_union_calls/no_reference"
readonly ref_path="${ref_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.gz"
readonly ref_type="vcf"

readonly phased_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
readonly phased_path="${phased_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly phased_type="vcf"

readonly out_dir="data/prephased/wes_union_calls/phase_conf"
readonly out_prefix="${out_dir}/ukb_shapeit5_whatshap_chr${chr}"
readonly out_type="mt"
readonly out="${out_prefix}.mt"

mkdir -p ${out_dir}



# combine parents and children in same vcf
# note: we can't combine data using bcftools beacuse some variants
# share the same position (but not same reference alleles)
module purge
set_up_hail
set_up_pythonpath_legacy
python3 ${hail_script} \
  --phased_path ${phased_path} \
  --phased_type ${phased_type} \
  --ref_path ${ref_path} \
  --ref_type ${ref_type} \
  --out_prefix ${out_prefix} \
  --out_type ${out_type}

