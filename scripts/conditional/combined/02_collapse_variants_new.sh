#!/usr/bin/env bash
#
# @description: amalgamate variants by phase to infer knockouts by genes.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=collapse_variants
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/collapse_variants.log
#SBATCH --error=logs/collapse_variants.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 4
#SBATCH --array=21
#
#
#$ -N collapse_variants
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/collapse_variants.log
#$ -e logs/collapse_variants.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 15
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/combined/02_collapse_variants.py"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/prefilter/pp90"
readonly out_dir="data/mt/dosages/pp90"
# in parameters
readonly in_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.mt"
readonly in_type="mt"
# prefix for indiviual genes and final merged file
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_max_ds"
readonly out_type="vcf"

readonly csqs_category="pLoF,damaging_missense"

mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
  --chrom ${chr} \
  --input_path ${in_prefix} \
  --input_type ${in_type} \
  --out_prefix ${out_prefix} \
  --out_type ${out_type} \
  --csqs_category ${csqs_category}

