#!/usr/bin/env bash
#
# Extract genetic coordinates for significant genes in primary analysis.
#
#$ -N gene_positions
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/gene_positions.log
#$ -e logs/gene_positions.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qf
#$ -t 1-80
#$ -tc 10

#set -o errexit
#set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/common/01_gene_positions.R"

readonly genes="data/genes/220310_ensgid_grch38_pos.tsv.gz"
readonly out_dir="data/conditional/common/intervals"
readonly pheno_dir="data/phenotypes"
readonly in_prefix="ukb_eur_wes_200k"
readonly maf="0to5e-2"


submit_binary_intervals()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "binary"
}

submit_cts_intervals()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "cts"
}

submit_intervals()
{
  set_up_rpy
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  local step2_dir="data/saige/output/${trait}/step2"
  
  local in_file="${step2_dir}/${in_prefix}_maf${maf}_${phenotype}_${annotation}.txt.gz"
  local out_prefix="${out_dir}/${in_prefix}_maf${maf}_${phenotype}_${annotation}"
  
  mkdir -p ${out_dir}
  
  set -x
  Rscript "${rscript}" \
    --in_spa_file "${in_file}" \
    --coordinates "${genes}" \
    --flanking_bp 0 \
    --fdr_cutoff "0.25" \
    --out_prefix "${out_prefix}"
  set +x
}

submit_binary_intervals "pLoF_damaging_missense"
submit_cts_intervals "pLoF_damaging_missense"
submit_binary_intervals "pLoF"
submit_cts_intervals "pLoF"
submit_binary_intervals "damaging_missense"
submit_cts_intervals "damaging_missense"
submit_binary_intervals "synonymous"
submit_cts_intervals "synonymous"




