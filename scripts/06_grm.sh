#!/usr/bin/env bash
#
# Create sparse genetic relatedness matrix using genotyped/imputed data.
#
#$ -N grm
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/grm.log
#$ -e logs/grm.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly grm_script="scripts/_create_grm.sh"
readonly mrg_script="scripts/_merge_grm.sh"
readonly fit_script="scripts/_fit_grm.sh"

readonly spark_dir="data/tmp/spark"
readonly out_dir="data/saige/grm/input/chunks"
readonly out_prefix="${out_dir}/ukb_wes_200k_sparse_autosomes"
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'


mkdir -p ${out_dir}

create_grm() {
  local out_grm_prefix=${1}
  local out_type="mt"
  set -x
  qsub -N "_create_grm" \
    -o "logs/_create_grm.log" \
    -e "logs/_create_grm.errors.log" \
    -t 1-22 \
    -q "short.qc" \
    -pe shmem 5 \
    "${grm_script}" \
    "${out_grm_prefix}" \
    "${out_type}" \
    "${final_sample_list}" 
  set +x
}

merge_grm() {
  local in_merge_prefix=${1}
  local out_merge_prefix=${2}
  local in_type="mt"
  local out_type="plink"
  set -x
  qsub -N "_merge_grm" \
    -o "logs/_merge_grm.log" \
    -e "logs/_merge_grm.errors.log" \
    -hold_jid "_create_grm" \
    -t 1 \
    -q "short.qc" \
    -pe shmem 4 \
    "${mrg_script}" \
    "${in_merge_prefix}" \
    "${in_type}"\
    "${out_merge_prefix}"\
    "${out_type}"
  set +x
}

fit_grm() {
  local plink_file=${1}
  local out_fit_prefix=${2}
  set -x
  qsub -N "_create_grm" \
    -o "logs/_create_grm.log" \
    -e "logs/_create_grm.errors.log" \
    -hold_jid "_merge_grm" \
    -q "short.qc" \
    -pe shmem 10 \
    "${grm_script}" \
    "${plink_file}" \
    "${out_fit_prefix}"
  set +x
}


create_grm ${out_prefix}








