#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=co_occurence_by_phenotype
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/co_occurence_by_phenotype.log
#SBATCH --error=logs/co_occurence_by_phenotype.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=16

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly rscript="scripts/post_hoc/co_occurence/02_co_occurence_by_phenotype.R"

# annotation for knockouts to be extractaed
readonly annotation="pLoF_damaging_missense"
# phenotype files to use (either tte or standard)
readonly in_dir="data/phenotypes"
readonly path_header="${in_dir}/dec22_phenotypes_binary_200k_header.tsv"
readonly path_phenotypes="${in_dir}/dec22_phenotypes_binary_200k.tsv.gz"
readonly path_tte_phenotypes="${in_dir}/tte_matrix_176k.txt.gz"

readonly out_dir="data/knockouts/alt/pp90/co_occurence6"
readonly out_prefix_tte="${out_dir}/co_occurence_by_tte_phenotype_chr${chr}"
readonly out_prefix="${out_dir}/co_occurence_by_phenotype_chr${chr}"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
  --chrom ${chr} \
  --annotation ${annotation} \
  --path_header ${path_header} \
  --out_prefix ${out_prefix} \
  --path_phenotypes ${path_phenotypes}
  #--path_phenotypes ${path_tte_phenotypes} \
  #--convert_tte_to_bool \


