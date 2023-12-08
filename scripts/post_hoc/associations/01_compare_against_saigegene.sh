#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=compare_against_saigegene
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/compare_against_saigegene.log
#SBATCH --error=logs/compare_against_saigegene.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

source utils/bash_utils.sh

readonly rscript="scripts/post_hoc/associations/01_compare_against_saigegene.R"
set_up_rpy && Rscript "${rscript}"

