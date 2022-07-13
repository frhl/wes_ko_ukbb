#!/usr/bin/env bash
#
#$ -N eval_chunks_by_method
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/eval_chunks_by_method.log
#$ -e logs/eval_chunks_by_method.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -V


source utils/bash_utils.sh

readonly rscript="scripts/phasing/07_eval_chunks_by_method.R"
readonly main_dir="data/phased/wes_union_calls/chunks"
readonly out_dir="derived/phased/validation"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_phasing"
readonly qc_sites="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
    --master_chunk_dir "${main_dir}" \
    --out_prefix "${out_prefix}" \
    --sites "${qc_sites}" \
    --img_height 5 \
    --img_width 8

