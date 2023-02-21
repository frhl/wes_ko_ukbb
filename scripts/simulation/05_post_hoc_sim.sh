#!/usr/bin/env bash
#
#$ -N post_hoc_sim
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/post_hoc_sim.log
#$ -e logs/post_hoc_sim.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/simulation/05_post_hoc_sim.R"

readonly in_dir="data/simulation/saige/step2/binary"
readonly out_dir="data/simulation/combined"
mkdir -p ${out_dir}

pattern_seed="seed10"
pattern_var="var_0.10"
out_prefix="${out_dir}/seed100_combined_sim"

module purge
set_up_rpy
Rscript ${rscript} \
  --input_dir "${in_dir}" \
  --out_prefix "${out_prefix}" \
  --seed_regex "${pattern_seed}" \
  --var_regex "${pattern_var}"


