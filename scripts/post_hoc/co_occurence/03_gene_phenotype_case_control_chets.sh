#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=gene_phenotype_case_control_chets
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/gene_phenotype_case_control_chets.log
#SBATCH --error=logs/gene_phenotype_case_control_chets.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 4


source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly rscript="scripts/post_hoc/03_gene_phenotype_case_control_chets.R"

readonly co_occurence_dir="data/knockouts/alt/pp90/co_occurence5"
readonly co_occurence_file="${co_occurence_dir}/co_occurence_by_phenotype_chrCHR.txt.gz"


readonly out_dir="data/knockouts/alt/pp90/co_occurence5"
readonly out_prefix="${out_dir}/co_occurence_collapsed_pLoF_damaging_missense"

mkdir -p ${out_dir}

## Note: remember to include chr1 in the R-script!
set_up_rpy
Rscript "${rscript}" \
  --co_occurence_file ${co_occurence_file} \
  --out_prefix ${out_prefix}


