#!/usr/bin/env bash

#SBATCH --account=lindgren.prj
#SBATCH --job-name=write_gene_ko
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/write_gene_ko.log
#SBATCH --error=logs/write_gene_ko.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --array=1-21


set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly rscript="scripts/survival/01_write_gene_ko.R"

# VCF used to get wildtype (hom ref) samples.
readonly vcf_dir="data/knockouts/alt/pp90/combined"
readonly vcf="${vcf_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense.vcf.bgz"



extract_genes_by_annotation() {
  local annotation=${1}
  local in_dir="data/knockouts/alt/pp90/combined"
  local in_file="${in_dir}/ukb_eur_wes_200k_chr${chr}_${annotation}_all.tsv.gz"
  local out_dir="data/survival/knockouts/${annotation}/wo_chrom_prefix"
  local out_prefix="${out_dir}/ukb_eur_wes_200k_${annotation}"
  local sample_file="${out_prefix}.samples.txt"
  if [ -f "${in_file}" ]; then
    mkdir -p ${out_dir}
    bcftools query -l ${vcf} > ${sample_file}
    Rscript ${rscript} \
      --in_file ${in_file} \
      --sample_file ${sample_file} \
      --out_prefix ${out_prefix}
    rm ${sample_file}
  else
    >&2 echo "${in_file} does not exist."
  fi

 }

set_up_rpy
module load BCFtools/1.12-GCC-10.3.0
extract_genes_by_annotation "pLoF_damaging_missense"
#extract_genes_by_annotation "pLoF"



