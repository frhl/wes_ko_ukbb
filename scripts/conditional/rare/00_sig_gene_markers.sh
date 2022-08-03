#!/usr/bin/env bash
#
#$ -N sig_gene_markers
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/sig_gene_markers.log
#$ -e logs/sig_gene_markers.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 1-80

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/rare/00_sig_gene_markers.R"

# parameters
readonly maf="0to5e-2"
readonly min_mac=4
readonly use_prs=1

# I/O
readonly genes="data/genes/220310_ensgid_grch38_pos.tsv.gz"
readonly out_dir="data/conditional/rare/combined/genes/min_mac${min_mac}"
readonly vep_dir="data/mt/vep"

readonly pheno_dir="data/phenotypes"
readonly in_prefix="ukb_eur_wes_200k"
readonly vep_path="${vep_dir}/ukb_eur_wes_200k_csqs_chrCHR.tsv.gz" # CHR to be gsubbed in Rscript

mkdir -p ${out_dir}

get_sig_gene_markers_binary()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_sig_genes "${annotation}" "${phenotype}" "binary"
}

get_sig_gene_markers_cts()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${SGE_TASK_ID}q;d" ${pheno_list} )
  submit_sig_genes "${annotation}" "${phenotype}" "cts"
}

submit_sig_genes()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  local step2_dir="data/saige/output/${trait}/step2/min_mac${min_mac}"
 
  # path to PRS files versus non PRS files are different 
  local in_file_loco="${step2_dir}/${in_prefix}_maf${maf}_${phenotype}_${annotation}_locoprs.txt.gz"
  local in_file_std="${step2_dir}/${in_prefix}_maf${maf}_${phenotype}_${annotation}.txt.gz"
  local out_prefix_std="${out_dir}/${in_prefix}_maf${maf}_${phenotype}_${annotation}"

  if [ "${use_prs}" -eq "1" ] & [ -f "${in_file_loco}" ]; then
    echo "Note: Using PRS results from ${phenotype}"
    local in_file=${in_file_loco}
  else
    echo "Warning: No PRS for ${phenotype}. Using standard results instead."
    local in_file=${in_file_std}
  fi

  local out_prefix=${out_prefix_std}
  set_up_rpy
  Rscript "${rscript}" \
    --vep_path_with_CHR "${vep_path}" \
    --in_spa_file "${in_file}" \
    --out_prefix "${out_prefix}" \
    --phenotype "${phenotype}"

}


get_sig_gene_markers_binary "pLoF_damaging_missense"




