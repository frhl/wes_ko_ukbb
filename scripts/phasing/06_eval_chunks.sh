#!/usr/bin/env bash
#
#$ -N eval_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/eval_chunks.log
#$ -e logs/eval_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -t 20
#$ -V


source utils/bash_utils.sh

readonly rscript="scripts/phasing/06_eval_chunks.R"

#readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly main_dir="data/phased/wes_union_calls/chunks/final"

END=22
in_dir=""
for chr in $(seq 1 $END); do 
  dir="${main_dir}/ukb_eur_wes_union_calls_200k_chr${chr}-16xshort.qe"
  in_dir="${in_dir},${dir}"
done

#readonly dir1="${main_dir}/ukb_wes_union_calls_200k_chr20-16xshort.qa"
#readonly dir2="${main_dir}/ukb_wes_union_calls_200k_chr20-24xshort.qa"
#readonly dir3="${main_dir}/ukb_wes_union_calls_200k_chr20-16xlong.qc"
#readonly in_dir="${dir1},${dir2},${dir3}"

readonly in_dir_merged="data/phased/wes_union_calls/merged"

readonly out_dir="derived/phased/wes_union_calls/qc"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls"

mkdir -p ${out_dir}


set_up_rpy
set -x
Rscript ${rscript} \
    --in_dir "${in_dir}" \
    --in_dir_merged "${in_dir_merged}" \
    --out_prefix "${out_prefix}" \
    --img_height 12
set +x

