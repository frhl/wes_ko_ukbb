#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_singleton_cat
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_singleton_cat.log
#SBATCH --error=logs/get_singleton_cat.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly rscript="scripts/post_hoc/variants/09_get_singleton_cat.R"

readonly vep_dir="data/vep/vep95/worst_csqs"
readonly vep_file="${vep_dir}/UKB.chrCHR.exome_array.variants_only.vep95.csqs.worst_csq_by_gene_canonical.original.txt.gz"

readonly in_dir="data/mt/prefilter/pp90"
readonly in_file="${in_dir}/ukb_wes_union_calls_200k_chrCHR.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.alt_alleles.txt.gz"

readonly out_dir="data/mt/singleton_category_count"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.from_mt.alt_alleles.singleton_category"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --vep_file ${vep_file} \
  --in_file ${in_file} \
  --out_prefix ${out_prefix}



