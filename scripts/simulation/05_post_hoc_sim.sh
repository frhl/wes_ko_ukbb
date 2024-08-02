#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=post_hoc_sim
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/post_hoc_sim.log
#SBATCH --error=logs/post_hoc_sim.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/simulation/05_post_hoc_sim.R"

readonly in_dir="data/simulation/saige/step2/binary"
readonly out_dir="data/simulation/combined"
mkdir -p ${out_dir}

#pattern_var="h2_0\\.00" 
# Simualate traits with no effect
pattern_seed="seed40"
pattern_var="K0.01_"
ac_allele2_cutoff=10
out_prefix="${out_dir}/seed40_combined_sim_null_ac10"

module purge
set_up_rpy
Rscript ${rscript} \
  --input_dir "${in_dir}" \
  --out_prefix "${out_prefix}" \
  --input_pattern "${pattern_seed}" \
  --ac_allele2_cutoff "${ac_allele2_cutoff}"



