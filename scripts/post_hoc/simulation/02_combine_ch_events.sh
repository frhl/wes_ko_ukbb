#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=combine_ch_events
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/combine_ch_events.log
#SBATCH --error=logs/combine_ch_events.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly array_idx=$( get_array_task_id )

readonly rscript="scripts/post_hoc/simulation/02_combine_ch_events.R"

readonly n_samples="176935"

readonly in_dir="data/simulation/sim_ch_events/2401"
readonly in_regex="ch_events"

readonly out_dir="data/simulation/sim_ch_events/2401"
readonly out_prefix="${out_dir}/combined_simulation_n${n_samples}"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --in_dir ${in_dir} \
  --in_regex ${in_regex} \
  --out_prefix ${out_prefix}



