#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=combine_vep_with_plink_frqx
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/combine_vep_with_plink_frqx.log
#SBATCH --error=logs/combine_vep_with_plink_frqx.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/variants/01_combine_vep_with_plink_frqx.R"

readonly ac_dir_before_pp_filter="data/mt/annotated/old"
readonly ac_path_before_pp_filter="${ac_dir_before_pp_filter}/ukb_wes_union_calls_200k_chrCHR.frqx"
readonly ac_dir_after_pp_filter="data/mt/prefilter/pp90"
readonly ac_path_after_pp_filter="${ac_dir_after_pp_filter}/ukb_wes_union_calls_200k_chrCHR.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.frqx"

readonly vep_dir="data/vep/vep95/worst_csqs/"
readonly vep_path="${vep_dir}/UKB.chrCHR.exome_array.variants_only.vep95.csqs.worst_csq_by_gene_canonical.original.txt.gz"

readonly out_dir="data/vep/counts"
readonly out_prefix="${out_dir}/UKB.exome_array.variants.vep95.worst_csq_by_gene_canonical.original.counts"


mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --AC_path_before_pp_filter ${ac_path_before_pp_filter} \
  --AC_path_after_pp_filter ${ac_path_after_pp_filter} \
  --vep_path ${vep_path} \
  --out_prefix ${out_prefix}









