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
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/04_merge_freqx_with_vep.R"
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

# VEP version
readonly vep_version="105"
readonly vep_dir="data/vep/vep${vep_version}/worst_csqs"
readonly vep_path="${vep_dir}/UKB.chr${chr}.exome_array.variants_only.vep${vep_version}.csqs.worst_csq_by_gene_canonical.original.txt.gz"

## prior to filtering
#readonly ac_dir="data/phased/wes_union_calls/200k/shapeit5/plink"
#readonly ac_path="${ac_dir}/ukb_wes_union_calls_200k_chr${chr}.frqx"

#readonly out_dir="data/mt/prefilter/pp90/summary/after_phasing/ve"
#readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/vep_freqx_summary/vep${vep_version}"
#readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.vep${vep_version}"

## after filtering
readonly ac_dir="data/mt/prefilter/pp90"
readonly ac_path="${ac_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.frqx"

readonly out_dir="data/mt/prefilter/pp90/summary"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.vep${vep_version}"


mkdir -p ${out_dir}


set_up_rpy
Rscript ${rscript} \
  --AC_path ${ac_path} \
  --vep_spliceAI_processed ${vep_path} \
  --force_without_spliceAI \
  --out ${out_prefix}




