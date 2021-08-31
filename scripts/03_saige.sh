#!/usr/bin/env bash
#
#
#$ -N saige_hail
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/saige_hail.log
#$ -e logs/saige_hail.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 22

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

SGE_TASK_ID=22

# directories
readonly in_dir_phased="data/phased"
readonly in_dir_unphased="data/unphased/unfiltered"
readonly vep_dir="data/vep/full"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/saige/input"

# hail script
readonly hail_script="utils/hail_plink_export.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly in_phased="${in_dir_phased}/ukb_wes_200k_phased_chr${chr}.1of1.vcf.gz"
readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"
readonly out="${out_prefix}.plink"

# setup hail
set_up_hail

SECONDS=0
#mkdir -p ${out_dir}
#python3 ${hail_script} \
#    --chrom ${chr} \
#    --input_path ${in_unphased} \
#    --input_type "mt" \
#    --get_europeans \
#    --missing 0.05 \
#    --out_prefix ${out_prefix} \
#    --out_type "plink"


Rscript createSparseGRM.R	\
	--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
	--nThreads=4  \
	--outputPrefix=./output/sparseGRM	\
	--numRandomMarkerforSparseKin=1000	\
	--relatednessCutoff=0.125


print_update "Finished running HAIL for chr${chr}" "${SECONDS}"

