#!/usr/bin/env bash

#SBATCH --account=lindgren.prj
#SBATCH --job-name=write_forest_dat_table
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/write_forest_dat_table.log
#SBATCH --error=logs/write_forest_dat_table.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly rscript="scripts/survival/03_write_forest_dat_table.R"

readonly in_dir="/well/lindgren-ukbb/projects/ukbb-11867/samvida/for_fred_ko_project/2304_analyses/results/ref_group_het"
readonly in_dir_cond="${in_dir}/unconditioned"
readonly in_dir_uncond="${in_dir}/unconditioned"

readonly out_dir="data/survival/tables"
readonly out_prefix="${out_dir}/230609_forest_table"
mkdir -p ${out_dir}

readonly n_tests=266560 # (952 * 280)

set_up_rpy
Rscript ${rscript} \
  --in_dir_cond ${in_dir_cond} \
  --in_dir_uncond ${in_dir_uncond} \
  --out_prefix ${out_prefix} \
  --n_tests ${n_tests}





