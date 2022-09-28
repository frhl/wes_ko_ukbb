#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_ld
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calc_ld.log
#SBATCH --error=logs/calc_ld.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 5
#SBATCH --requeue


source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/04_calc_ld.R"

readonly bed_dir="data/prs/hapmap/ld/unrel_kin_eur_10k"
readonly bed_file="${bed_dir}/short_merged_ukb_hapmap_rand_10k_eur.bed"

readonly out_dir="data/prs/hapmap/ld/matrix_unrel_kin"
readonly out_prefix="${out_dir}/ld_matrix"

mkdir -p ${out_dir}

export OPENBLAS_NUM_THREADS=1 # avoid two levels of parallelization

fit_ld_matrix() 
{
  SECONDS=0
  set_up_rpy
  set -x
  Rscript "${rscript}" \
   --bed "${bed_file}" \
   --out_prefix "${out_prefix}"
  set +x
  log_runtime ${SECONDS}
}


fit_ld_matrix



