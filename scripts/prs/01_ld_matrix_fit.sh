#!/usr/bin/env bash
#
#$ -N ld_matrix_fit
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ld_matrix_fit.log
#$ -e logs/ld_matrix_fit.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q long.qc@@long.hge
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly rscript="scripts/prs/00_bed_gen.py"
readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly out_dir="data/prs/hapmap"
readonly out_prefix="${out_dir}/long_ukb_ld_hapmap_rand_10k_eur"

mkdir -p ${out_dir}

if [ ! -f "${out_prefix}" ]; then
  set_up_rpy
  set -x
  Rscript "${rscript}" \
     --chrom "AUTOSOMES" \
     --dataset "imp" \
  set +x
else
  print_update "file ${out} already exists. Skipping!"
fi




