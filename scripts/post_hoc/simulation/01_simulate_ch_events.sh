#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=simulate_ch_events
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/simulate_ch_events.log
#SBATCH --error=logs/simulate_ch_events.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=2-10

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly array_idx=$( get_array_task_id )

readonly rscript="scripts/post_hoc/simulation/01_simulate_ch_events.R"

readonly n_samples="176935"
readonly seed="${array_idx}"

readonly out_dir="data/simulation/sim_ch_events/2401"
readonly out_prefix="${out_dir}/ch_events_seed${seed}_n${n_samples}"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --out_prefix ${out_prefix} \
  --n_samples ${n_samples} \
  --seed ${seed}



