#!/usr/bin/env bash
#
#$ -N export_ko_geneset
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/export_ko_geneset.log
#$ -e logs/export_ko_geneset.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qf

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly rscript="scripts/geneset/03_export_ko_geneset.R"

readonly in_dir="data/geneset/knockouts"
readonly in_vcf="${in_dir}/ukb_eur_wes_200k_maf0to5e-2_pLoF_damaging_missense_combined.vcf.bgz"

readonly out_dir="data/geneset/knockouts"
readonly out="${out_dir}/ukb_eur_wes_200k_maf0to5e-2_pLoF_damaging_missense_msigdb_h.saige"

readonly bridge="/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz"

mkdir -p ${out_dir}

# Generate SAIGE-GENE+ Group file consequence
# annotations (SAIGE version > 0.99.2)
module purge
set_up_rpy
Rscript ${rscript} \
  --in_vcf "${in_vcf}" \
  --output_path "${out}" \
  --bridge "${bridge}"


