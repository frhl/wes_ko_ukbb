#!/usr/bin/env bash
#
# Annotate variants
#
# Author: Frederik Lassen (2021-06-25)
#
#$ -N vep
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/vep.log
#$ -e logs/vep.log
#$ -P lindgren.prjc
#$ -pe shmem 8
#$ -q short.qc
#$ -t 22
#$ -V

set -o errexit
set -o nounset
module purge 

# Set variables
readonly chr=${SGE_TASK_ID}
readonly RAW_ROOT="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/phased"
readonly RAW_FILE="/ukb_wes_200k_phased_chr${chr}.1of1.vcf.gz"
readonly TMP_FILE="/ukb_wes_200k_phased_tmp_chr${chr}.1of1.vcf.gz"

readonly OUT_ROOT="derived/vep/output/test"
readonly OUT_FILE1="/ukb_wes_200k_vep_chr${chr}.vcf" # direct output of VEP

# extract variant information
module load BCFtools/1.9-foss-2018b # for extracting variant information

bcftools query -f '%CHROM %POS %ID %REF %ALT %AC %QUAL\n' "${RAW_ROOT}${RAW_FILE}" -o "${OUT_ROOT}${TMP_FILE}"

# Load modules (load order is important)
module purge
module load EnsEMBLCoreAPI/96.0-r20190601-foss-2019a-Perl-5.28.1 # required for LOFTEE
module load VEP/95.0-foss-2018b-Perl-5.28.0 # required FOR VEP
module load samtools/1.8-gcc5.4.0 # required for LOFTEE

export PERL5LIB=$PERL5LIB:/well/lindgren/flassen/software/VEP/plugins_grch38/

## run VEP 95 on temporary file
vep --input_file "${OUT_ROOT}${TMP_FILE}" \
--dir_cache /well/lindgren/flassen/software/VEP/vep95/GRCh38 \
--assembly GRCh38 \
--cache \
--species "homo_sapiens" \
--symbol \
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

# change names
mv "${OUT_ROOT}${OUT_FILE1}.tmp" "${OUT_ROOT}${OUT_FILE1}"
chmod -w "${OUT_ROOT}${OUT_FILE1}" 

