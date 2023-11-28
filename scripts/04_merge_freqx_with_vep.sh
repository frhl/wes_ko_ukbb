#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=merge_freqx_with_vep
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/merge_freqx_with_vep.log
#SBATCH --error=logs/merge_freqx_with_vep.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=22

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/04_merge_freqx_with_vep.R"
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly ac_dir="data/mt/prefilter/pp90"
readonly ac_path="${ac_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.frqx"

readonly vep_version="105"
readonly vep_dir="data/vep/vep${vep_version}/worst_csqs"
readonly vep_path="${vep_dir}/UKB.chr22.exome_array.variants_only.vep${vep_version}.csqs.worst_csq_by_gene_canonical.original.txt.gz"

readonly out_dir="data/mt/prefitler/pp90"
readonly out_prefix="${ac_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.vep${vep_version}"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --AC_path ${ac_path} \
  --vep_spliceAI_processed ${vep_path} \
  --force_without_spliceAI \
  --out ${out_prefix}




