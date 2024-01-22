#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=count_plotting
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/count_plotting.log
#SBATCH --error=logs/count_plotting.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/04_count_plotting.R"

# VEP version
readonly vep_version="95"

readonly in_dir_before_pp_cutoff="data/phased/wes_union_calls/200k/shapeit5/vep_freqx_summary/vep${vep_version}"
readonly in_dir_after_pp_cutoff="data/phased/wes_union_calls/200k/shapeit5/filter_by_pp"
readonly in_dir_after_pp_cutoff_with_qc="data/mt/prefilter/pp90/summary"
readonly out_dir="data/mt/prefilter/pp90/summary"
readonly out="${out_dir}/ukb_wes_union_calls_200k.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.vep${vep_version}.pre_post_filter.pdf"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --count_directory ${in_dir_after_pp_cutoff_with_qc} \
  --count_directory_other ${in_dir_after_pp_cutoff} \
  --count_directory_before_pp_cutoff ${in_dir_before_pp_cutoff} \
  --out ${out}




