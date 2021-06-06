#!/bin/bash

#$ -cwd
#$ -N vep10
#$ -o vep10.log
#$ -e vep10.errors.log
#$ -q short.qf
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -t 10

# Set variables
CHR=${SGE_TASK_ID}
RAW_ROOT=/well/lindgren/UKBIOBANK/nbaya/resources/ukb_wes_200k_inliers_split_filtered/
RAW_FILE=ukb_wes_200k_inliers_split_filtered_hail_chr${CHR}.vcf.bgz
TMP_FILE=ukb_wes_200k_inliers_split_filtered_chr${CHR}.vcf

OUT_ROOT=/well/lindgren/UKBIOBANK/flassen/projects/KO/WES200K/derived/vep/
OUT_FILE1=ukb_wes_200k_inliers_split_filtered_chr${CHR}.vcf # direct output of VEP
OUT_FILE2=ukb_wes_200k_vep_chr${CHR}.txt # output of merging RAW_FILE with OUT_FILE1

VEP=/well/lindgren/users/uyi901/software/ensembl-vep/vep

# Load modules
module load BCFtools/1.9-foss-2018b # for extracting variant information
module load ensembl-tools/94 # module required for VEP
module load BioPerl/1.7.2-foss-2018b-Perl-5.28.0 # required for LOFTEE
module load Bio-DB-HTS/2.11-foss-2018b-Perl-5.28.0 # required for LOFTEE
module load samtools/1.8-gcc5.4.0 # required for LOFTEE

# temporary file that only includes SNP level information
#bcftools query -f '%CHROM  %POS %ID  %REF  %ALT %AF %QUAL\n' "${RAW_ROOT}${RAW_FILE}" -o "${OUT_ROOT}${TMP_FILE}"
bcftools query -f '%CHROM  %POS %ID  %REF  %ALT\n' "${RAW_ROOT}${RAW_FILE}" -o "${OUT_ROOT}${TMP_FILE}"

export PERL5LIB=$PERL5LIB:/well/lindgren/flassen/software/VEP/plugins_grch38
export PERL5LIB=$PERL5LIB:/well/lindgren/flassen/software/VEP/vep_plugins/loftee

###########
# Run VEP #
###########

## run VEP 94 on temporary file
$VEP --input_file "${OUT_ROOT}${TMP_FILE}" \
--dir_cache /well/lindgren/flassen/software/VEP/GRCh38 \
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
--dir_plugins /well/lindgren/flassen/software/VEP/plugins_grch38 \
--plugin dbNSFP,/well/lindgren/flassen/software/VEP/dbNSFP/GRCh38/dbNSFP4.1a_grch38.gz,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred \
--fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CDS_position,Protein_position,EXON,INTRON,LoF_flags,LoF_filter,LoF,SIFT,PolyPhen,SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred" \
--output_file "$OUT_ROOT"TMPXX_VEP"$CHR" \
--force_overwrite \
--no_stats \
--verbose \
--offline

mv "$OUT_ROOT"TMP_VEP"$CHR" "${OUT_ROOT}${OUT_FILE1}"





