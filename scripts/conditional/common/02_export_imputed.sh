#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_imputed
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/export_imputed.log
#SBATCH --error=logs/export_imputed.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=6
#
#
#$ -N export_imputed
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/export_imputed.log
#$ -e logs/export_imputed.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -t 22
#$ -V


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly hail_script="scripts/conditional/common/02_export_imputed.py"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

readonly min_maf=0.01
readonly min_info=0.8
readonly missing=0.10

readonly out_dir="data/unphased/imputed/common"
readonly out_prefix="${out_dir}/ukb_imp_200k_common"
readonly out_type="mt"

mkdir -p ${out_dir}

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
     --chrom ${chr} \
     --min_maf ${min_maf} \
     --min_info ${min_info} \
     --missing ${missing} \
     --extract ${final_sample_list} \
     --out_prefix ${out_prefix} \
     --out_type ${out_type} \
     && print_update "Finished filtering imputed genotypes ${out_prefix}" ${SECONDS} \
     || raise_error "Filtering imputed genotypes for for ${out_prefix} failed!"



