#!/usr/bin/env bash
#
# @description evaluate quality of chunks by methods
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=eval_chunks_by_method
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/eval_chunks_by_method.log
#SBATCH --error=logs/eval_chunks_by_method.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

source utils/bash_utils.sh

readonly rscript="scripts/phasing/07_eval_chunks_by_method.R"
readonly main_dir="data/phased/wes_union_calls/chunks"
readonly out_dir="data/phased/validation"
readonly out_prefix="${out_dir}/220713_ukb_eur_wes_union_calls_phasing"
readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
    --master_chunk_dir "${main_dir}" \
    --out_prefix "${out_prefix}" \
    --sites "${wes_variants}" \
    --img_height 5 \
    --img_width 8

