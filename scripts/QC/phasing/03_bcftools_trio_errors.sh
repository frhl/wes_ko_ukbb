#!/usr/bin/env bash
#
#$ -N trio_switch_rate
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/trio_switch_rate.log
#$ -e logs/trio_switch_rate.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 21

source utils/bash_utils.sh
source utils/hail_utils.sh
module load BCFtools/1.9-foss-2018b

readonly in_dir_phased="data/phased/test-phasing"
readonly in_dir_unphased="data/unphased/post-qc"
readonly out_dir="data/mendel"
readonly fam_dir="/well/lindgren/UKBIOBANK/nbaya/resources"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly in_phased="${in_dir_phased}/ukb_wes_phased_non_singleton_chr${chr}-24xlong.qc-v4.2.2.vcf.gz"
#readonly in_phased="${in_dir_phased}/ukb_wes_200k_refphased_chr${chr}.vcf.gz"
#readonly in_phased="${in_dir_phased}/ukb_wes_200k_all_phased_chr${chr}.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly fam_file="${fam_dir}/ukb11867_pedigree.fam"


#readonly out_prefix="${out_dir}/211212_ukb_wes_trios_all_phased_chr${chr}"

export BCFTOOLS_PLUGINS="/apps/eb/skylake/software/BCFtools/1.9-foss-2018b/libexec/bcftools"

bcftools +trio-swich-rate "${in_phased}" -- -p "${fam_file}"

# calculate switch errors
#vcftools \
#  --gzvcf "${out_prefix}_estimation.vcf" \
#  --diff "${out_prefix}_transmission.vcf" \
#  --diff-switch-error \
#  --out "${out_prefix}"


# Calculate binominal confidence intervals
#set_up_rpy
#Rscript "${rscript}" \
#    --in_file "${out_prefix}.tsv.gz"\
#    --out_file "${out_prefix}_conf.tsv.gz"

