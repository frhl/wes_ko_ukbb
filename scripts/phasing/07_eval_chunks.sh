#!/usr/bin/env bash
#
#$ -N eval_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/eval_chunks.log
#$ -e logs/eval_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -V


source utils/bash_utils.sh

readonly rscript="scripts/phasing/07_eval_chunks.R"
readonly main_dir="data/phased/wes_union_calls/chunks"
readonly out_dir="data/phased/validation"
readonly out_prefix="${out_dir}/220713_ukb_eur_wes_union_calls"

readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
    --master_chunk_dir "${main_dir}" \
    --out_prefix "${out_prefix}" \
    --sites "${wes_variants}"\
    --img_height 12 \
    --img_width 8

