#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=simulate_ch_events
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/simulate_ch_events.log
#SBATCH --error=logs/simulate_ch_events.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --time=12:00:00

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly rscript="scripts/post_hoc/simulation/01_simulate_ch_events.R"


mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript}



