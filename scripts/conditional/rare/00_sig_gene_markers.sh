#!/usr/bin/env bash
#
# @description Takes gene hits that are significant in primary analysis
# and digs out all variants in the gene as annotated by VEP.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=sig_gene_markers
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/sig_gene_markers.log
#SBATCH --error=logs/sig_gene_markers.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 1
#SBATCH --array=1-320
#SBATCH --requeue


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh

readonly rscript="scripts/conditional/rare/00_sig_gene_markers.R"

# parameters
readonly min_mac=4
# Note: we generally don't want to use PRS here, since this is just a filter
# for how many genes to interrogate downstream. We will show PRS P-values regardless, 
# but would be nice to see how variants behave after conditionionon common/rare stuff.
readonly use_prs=0

# I/O
readonly genes="data/genes/220310_ensgid_grch38_pos.tsv.gz"
readonly out_dir="data/conditional/rare/combined/genes/2024/min_mac${min_mac}"
readonly vep_dir="data/mt/vep/worst_csq_by_gene_canonical"

readonly pheno_dir="data/phenotypes"
readonly in_prefix="ukb_eur_wes_200k"
readonly vep_path="${vep_dir}/ukb_eur_wes_union_calls_200k_chrCHR.tsv.gz" # CHR to be gsubbed in Rscript


# nominal significance P-cutoff
readonly genes_tested=952
readonly phenos_tested=311
readonly p_cutoff="$(python -c "print(0.05/(${genes_tested}))")"


readonly index=${SLURM_ARRAY_TASK_ID}

mkdir -p ${out_dir}

get_sig_gene_markers_binary()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/spiros_brava_phenotypes_binary_200k_header.tsv"
  local phenotype=$( sed "${index}q;d" ${pheno_list} )
  submit_sig_genes "${annotation}" "${phenotype}" "binary"
}

get_sig_gene_markers_cts()
{
  local annotation="${1?Error: Missing arg1 (annotation)}"
  local pheno_list="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  local phenotype=$( sed "${index}q;d" ${pheno_list} )
  submit_sig_genes "${annotation}" "${phenotype}" "cts"
}

submit_sig_genes()
{
  local annotation=${1?Error: Missing arg1 (consequence)}
  local phenotype=${2?Error: Missing arg2 (phenotype)}
  local trait=${3?Error: Missing arg3 (trait)}
  local step2_dir="data/saige/output/${trait}/step2/min_mac${min_mac}"

  if [ ! -z ${phenotype} ]; then
    if [ -f ${vep_path/CHR/21} ]; then

      # the possible paths that we can use
      local in_file_loco="${step2_dir}/${in_prefix}_${phenotype}_${annotation}_locoprs.txt.gz"
      local in_file_std="${step2_dir}/${in_prefix}_${phenotype}_${annotation}.txt.gz"
      local out_prefix_std="${out_dir}/${in_prefix}_${phenotype}_${annotation}"

      # depending on whether PRS is available select path
      if [ "${use_prs}" -eq "1" ] && [ -f "${in_file_loco}" ]; then
        echo "Note: Using PRS results from ${phenotype}"
        local in_file=${in_file_loco}
        local out_prefix="${out_prefix_std}_locoprs"
      else
        echo "Note: No PRS for ${phenotype}."
        local in_file=${in_file_std}
        local out_prefix="${out_prefix_std}"
      fi

      # Run R to get variants within genes
      set_up_rpy
      Rscript "${rscript}" \
        --vep_path_with_CHR "${vep_path}" \
        --in_spa_file "${in_file}" \
        --out_prefix "${out_prefix}" \
        --phenotype "${phenotype}" \
        --p_cutoff "${p_cutoff}"
    else 
      >&2 echo "${vep_path} (path) does not contain any files!."
    fi
 fi

}


get_sig_gene_markers_binary "pLoF_damaging_missense"




