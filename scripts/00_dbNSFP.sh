#!/usr/bin/env bash
#
# run vep to get dbNSFP annotations (currently not available through hail VEP)
#
#$ -N dbNSFP
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/dbNSFP.log
#$ -e logs/dbNSFP.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qe
#$ -t 1-21
#$ -V

set -o errexit
set -o nounset

source utils/hail_utils.sh

readonly chr=${SGE_TASK_ID}
readonly RAW_ROOT="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/unphased/post-qc"
readonly RAW_FILE="/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
readonly TMP_FILE="/ukb_wes_200k_filtered_tmp_chr${chr}.vcf.gz"

readonly OUT_ROOT="data/vep/full" #output"
readonly OUT_FILE1="/ukb_wes_200k_full_vep_chr${chr}.vcf" 

# Check if outfile already exists
if [ ! -f "${OUT_ROOT}${OUT_FILE1}" ]; then

  # extract variants in format that can be read by VEP
  if [ ! -f "${OUT_ROOT}${TMP_FILE}" ]; then
    module load BCFtools/1.9-foss-2018b
    bcftools query -f '%CHROM %POS %ID %REF %ALT %AC %QUAL\n' "${RAW_ROOT}${RAW_FILE}" -o "${OUT_ROOT}${TMP_FILE}"
  fi

  # Load modules (NOTE: load order is important)
  module purge
  set_up_vep

  ## run VEP 95 on temporary file
  vep --input_file "${OUT_ROOT}${TMP_FILE}" \
  --dir_cache /well/lindgren/flassen/software/VEP/vep95/GRCh38 \
  --assembly GRCh38 \
  --cache \
  --species "homo_sapiens" \
  --biotype \
  --canonical \
  --vcf \
  --pick \
  --total_length \
  --sift b \
  --polyphen b \
  --dir_plugins /well/lindgren/flassen/software/VEP/plugins_grch38/ \
  --plugin LoF,loftee_path:/well/lindgren/flassen/software/VEP/plugins_grch38/,human_ancestor_fa:/well/lindgren/flassen/software/VEP/loftee_human_ancestor/GRCh38/human_ancestor.fa.gz,conservation_file:/well/lindgren/flassen/software/VEP/loftee_human_ancestor/GRCh38/loftee.sql,gerp_bigwig:/well/lindgren/flassen/software/VEP/loftee_human_ancestor/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
  --plugin dbNSFP,/well/lindgren/flassen/software/VEP/dbNSFP/GRCh38/dbNSFP4.1a_grch38.gz,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationAssessor_score,MutationAssessor_pred,MutationTaster_score,MutationTaster_pred,PROVEAN_score,PROVEAN_pred,REVEL_score,CADD_raw,CADD_phred \
  --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CDS_position,Protein_position,EXON,INTRON,LoF_flags,LoF_filter,LoF,SIFT,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationAssessor_score,MutationAssessor_pred,MutationTaster_score,MutationTaster_pred,PROVEAN_score,PROVEAN_pred,REVEL_score,CADD_raw,CADD_phred" \
  --output_file "${OUT_ROOT}${OUT_FILE1}.tmp" \
  --force_overwrite \
  --no_stats \
  --verbose \
  --offline

  # change names and remove files
  mv "${OUT_ROOT}${OUT_FILE1}.tmp" "${OUT_ROOT}${OUT_FILE1}"
  mv "${OUT_ROOT}${OUT_FILE1}.tmp_warnings.txt" "${OUT_ROOT}/warnings.ukb_wes_200k_vep${chr}.vcf" 
  rm "${OUT_ROOT}${TMP_FILE}"
  chmod -w "${OUT_ROOT}${OUT_FILE1}"
else
  echo "${OUT_ROOT}${OUT_FILE1} already exists. Skipping... "
fi  


