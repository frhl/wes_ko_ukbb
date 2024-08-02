#!/usr/bin/env bash
#
# Extract genetic coordinates for significant genes in primary analysis.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=gene_positions
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/gene_positions.log
#SBATCH --error=logs/gene_positions.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-340

#$ -N gene_positions
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/gene_positions.log
#$ -e logs/gene_positions.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qa
#$ -t 1-330
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly task_id=$( get_array_task_id )

readonly rscript="scripts/conditional/common/01_gene_positions.R"
readonly rscript_prs="scripts/_check_prs_ok.R"

readonly min_mac=4

readonly genes="data/genes/220310_ensgid_grch38_pos.tsv.gz"
readonly out_dir="data/conditional/common/gene_positions/2402/min_mac${min_mac}"
readonly pheno_dir="data/phenotypes"
readonly in_prefix="ukb_eur_wes_200k"

readonly phenos_tested="311"
readonly genes_tested=952
readonly bonf_p_cutoff="$(python -c "print(0.05/(${genes_tested}*${phenos_tested}))")"
readonly nom_p_cutoff="$(python -c "print(0.05/(${genes_tested}))")"
echo "Using cutoff 'P < ${nom_p_cutoff}'."

mkdir -p ${out_dir}

submit_binary_intervals()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/dec22_phenotypes_binary_200k_header.tsv"
  local phenotype=$( sed "${task_id}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "binary"
}

submit_cts_intervals()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${task_id}q;d" ${pheno_list} )
  submit_intervals "${annotation}" "${phenotype}" "cts"
}

submit_intervals()
{
  set_up_rpy
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  local step2_dir="data/saige/output/${trait}/step2/min_mac${min_mac}"
  
  # Only use PRS if enabled and PRS file present
  local in_file_loco="${step2_dir}/${in_prefix}_${phenotype}_${annotation}_locoprs.txt.gz"
  local in_file_std="${step2_dir}/${in_prefix}_${phenotype}_${annotation}.txt.gz"
  #local out_prefix_loco="${out_dir}/${in_prefix}_${phenotype}_${annotation}_locoprs"
  local out_prefix_std="${out_dir}/${in_prefix}_${phenotype}_${annotation}"
  local prs_ok=$(Rscript ${rscript_prs} --phenotype ${phenotype})
  if [ "${use_prs}" -eq "1" ] & [ -f "${in_file_loco}" ] & [ "${prs_ok}" -eq "1" ]; then
    echo "Note: Using PRS results from ${phenotype}"
    local in_file=${in_file_loco}
  else
    echo "Warning: No PRS for ${phenotype}. Using standard results instead."
    local in_file=${in_file_std}
  fi
  local out_prefix=${out_prefix_std}

  if [ -f ${in_file} ]; then
    Rscript "${rscript}" \
      --in_spa_file "${in_file}" \
      --coordinates "${genes}" \
      --flanking_bp 0 \
      --out_prefix "${out_prefix}" \
      --phenotype "${phenotype}" \
      --p_cutoff "${nom_p_cutoff}"
  else
    >&2 echo "${in_file} does not exist. Skipping!"
  fi
}

# Use PRS when available
readonly use_prs=1

submit_binary_intervals "pLoF_damaging_missense"
#submit_cts_intervals "pLoF_damaging_missense"
#submit_binary_intervals "pLoF"
#submit_cts_intervals "pLoF"
#submit_binary_intervals "damaging_missense"
#submit_cts_intervals "damaging_missense"
#submit_binary_intervals "synonymous"
#submit_cts_intervals "synonymous"




