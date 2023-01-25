#!/usr/bin/env bash
#
# @description: amalgamate variants by phase to infer knockouts by genes.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=clinvar_alleles
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/clinvar_alleles.log
#SBATCH --error=logs/clinvar_alleles.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-19
#
#
#$ -N clinvar_alleles
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/clinvar_alleles.log
#$ -e logs/clinvar_alleles.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 1-22
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/06_clinvar_alleles.py"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/prefilter/pp90"
readonly out_dir="data/knockouts/alt/pp90/clinvar_alleles"
readonly in_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.mt"
readonly in_type="mt"

readonly clinvar_dir="/well/lindgren/flassen/ressources/clinvar/ftp"
readonly clinvar_path="${clinvar_dir}/clinvar_20230107.vcf.gz"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k_clinvar_chr${chr}"

mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --input_path ${in_prefix}\
   --input_type "mt" \
   --clinvar_path ${clinvar_path} \
   --out_prefix ${out_prefix}





