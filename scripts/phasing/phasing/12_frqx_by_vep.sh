#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=frqx_by_vep
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/frqx_by_vep.log
#SBATCH --error=logs/frqx_by_vep.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-21

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/phasing/phasing/12_frqx_by_vep.R"
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly vep_version="95"
readonly vep_dir="data/vep/vep${vep_version}/worst_csqs"
readonly vep_path="${vep_dir}/UKB.chr${chr}.exome_array.variants_only.vep${vep_version}.csqs.worst_csq_by_gene_canonical.original.txt.gz"

readonly ac_dir="data/phased/wes_union_calls/200k/shapeit5/filter_by_pp"
readonly ac_path="${ac_dir}/ukb_wes_union_calls_200k_chr${chr}.clean.frqx.gz"

readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/filter_by_pp"
readonly out_prefix="${ac_dir}/ukb_wes_union_calls_200k_chr${chr}.clean.frqx.vep${vep_version}"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --AC_path ${ac_path} \
  --vep_spliceAI_processed ${vep_path} \
  --force_without_spliceAI \
  --out ${out_prefix}

