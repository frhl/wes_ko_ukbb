#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=agreement
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/agreement.log
#SBATCH --error=logs/agreement.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#
#
#$ -N agreement
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/agreement.log
#$ -e logs/agreement.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

source utils/bash_utils.sh

readonly rscript="scripts/phasing/prephasing/07_agreement.R"
readonly input_path="data/phased/wes_union_calls/200k/calibration/ukb_shapeit5_whatshap_variants_chr21.PS.txt.gz"

readonly out_dir="data/prephased/wes_union_calls/"
readonly out_prefix="${out_dir}/221124_whatshap_s5_agreement_n10000.txt"

readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --input_path "${input_path}" \
  --n_samples 10000 \
  --seed 11415 \
  --sites "${wes_variants}" \
  --output_path "${out_prefix}" 






