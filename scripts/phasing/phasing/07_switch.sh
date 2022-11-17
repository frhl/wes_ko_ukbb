#!/usr/bin/env bash
#
# @description evaluate phasing quality stratified by MAF
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_ser
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/calc_ser.log
#SBATCH --error=logs/calc_ser.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2

source utils/bash_utils.sh

readonly rscript="scripts/phasing/11_calc_ser.R"
readonly ligated_dir="data/phased/wes_union_calls/200k/shapeit5/parents"
readonly out_dir="data/phased/validation"

readonly out_prefix="${out_dir}/221116_switch_error_rates"

readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --ligated_dir "${ligated_dir}" \
  --sites "${wes_variants}" \
  --out_prefix "${out_prefix}" 






