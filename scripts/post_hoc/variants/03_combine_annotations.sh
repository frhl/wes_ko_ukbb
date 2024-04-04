#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=combine_annotations
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/combine_annotations.log
#SBATCH --error=logs/combine_annotations.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 10
#SBATCH --requeue

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/post_hoc/03_combine_annotations.R"

readonly out_dir="data/knockouts/tables"
readonly out_prefix="${out_dir}/combined_annotations_by_sample.new"

mkdir -p ${out_dir}

set_up_rpy
Rscript "${rscript}" \
   --out_prefix "${out_prefix}"


