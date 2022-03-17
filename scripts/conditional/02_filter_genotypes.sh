#!/usr/bin/env bash
#
# extract genotypes in regions near genes that are significant in primary analysis
#
#$ -N filter_genotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/filter_genotypes.log
#$ -e logs/filter_genotypes.errors.log
#$ -P lindgren.prjc
#$ -q test.qc
#$ -t 1-40
#$ -V


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly bash_script="scripts/conditional/_filter_genotypes.sh"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

readonly padding=10e+6
readonly min_maf=0.01
readonly min_info=0.8

readonly in_dir="data/conditional/common/intervals"
readonly out_dir="data/conditional/common/intervals"
readonly pheno_dir="data/phenotypes"
readonly in_prefix="ukb_eur_wes_200k"
readonly maf="0to5e-2"

submit_binary_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "binary"
}

submit_cts_analysis()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "cts"
}

submit_intervals()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}

  local genes="${in_dir}/${in_prefix}_maf${maf}_${phenotype}_${annotation}.tsv.gz"
  local out_prefix="${out_dir}/${in_prefix}_maf${maf}_${phenotype}_${annotation}"

  mkdir -p ${out_dir}
  set -x
  qsub -N "_filter_genotypes_${1}" \
    -q "short.qc@@short.hge" \
    -t "${SGE_TASK_ID}" \
    -pe shmem 4 \
    "${bash_script}" \
    "${genes}" \
    "${final_sample_list}" \
    "${out_prefix}" \
    "${padding}" \
    "${min_maf}" \
    "${min_info}"
  set +x
}

submit_binary_analysis "pLoF_damaging_missense"
submit_cts_analysis "pLoF_damaging_missense"




