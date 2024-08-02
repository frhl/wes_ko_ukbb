#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_freqx_maf
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calc_freqx_maf.log
#SBATCH --error=logs/calc_freqx_maf.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/04_calc_freqx_maf.R"
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly vep_version="95"
readonly vep_dir="data/vep/vep${vep_version}/worst_csqs"
readonly vep_path="${vep_dir}/UKB.chr${chr}.exome_array.variants_only.vep${vep_version}.csqs.worst_csq_by_gene_canonical.original.txt.gz"

readonly ac_dir="data/mt/prefilter/pp90"
readonly ac_path="${ac_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.frqx"

readonly out_dir="data/mt/prefilter/pp90"
readonly out="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.frqx.maf.gz"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --vep_path ${vep_path} \
  --AC_path ${ac_path} \
  --out ${out}




