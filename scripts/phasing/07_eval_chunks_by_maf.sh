#!/usr/bin/env bash
#
#$ -N eval_chunks_by_maf
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/eval_chunks_by_maf.log
#$ -e logs/eval_chunks_by_maf.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V


source utils/bash_utils.sh

readonly rscript="scripts/phasing/07_eval_chunks_by_maf.R"
readonly ligated_dir="data/phased/wes_union_calls/ligated"
readonly out_dir="data/phased/validation"

readonly out_prefix="${out_dir}/220713_ligated"

readonly wes_variants="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
    --ligated_dir "${ligated_dir}" \
    --sites "${wes_variants}" \
    --out_prefix "${out_prefix}" 




