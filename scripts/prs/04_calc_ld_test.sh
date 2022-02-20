#!/usr/bin/env bash
#
#$ -N calc_ld_test
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calc_ld_test.log
#$ -e logs/calc_ld_test.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qc@@short.hga
#$ -V

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/prs/04_calc_ld_test.R"

readonly bed_dir="data/prs/hapmap/ld/unrel_eur_10k"
readonly bed_file="${bed_dir}/short_merged_ukb_hapmap_rand_10k_eur.bed"

readonly out_dir="data/prs/hapmap/ld/matrix"
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


fit_ld_matrix ${phenotype_binary}



