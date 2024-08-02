#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_collapsed_markers
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_collapsed_markers.log
#SBATCH --error=logs/get_collapsed_markers.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/combined/02_get_collapsed_markers.R"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly vep_dir="data/mt/prefilter/final_90_loftee"
readonly vep="${vep_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.csqs.txt.gz"

readonly rare_dir="data/conditional/rare/combined/mt"
readonly rare_file="${rare_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_markers.txt.gz"

readonly vcf_dir="data/mt/dosages_urv/pp90"
readonly vcf="${vcf_dir}/ukb_eur_wes_200k_chr${chr}_max_ds.vcf.bgz"

readonly out_dir="data/mt/dosages_urv/pp90"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_max_ds"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --path_worst_csq_by_gene_canonical ${vep} \
  --path_markers ${rare_file} \
  --path_vcf ${vcf} \
  --chrom ${chr} \
  --out_prefix ${out_prefix}






