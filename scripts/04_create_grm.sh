#!/usr/bin/env bash
#
#
#$ -N make_grm
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/make_grm.log
#$ -e logs/make_grm.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q short.qe

#set -o errexit
#set -o nounset

module purge
source utils/bash_utils.sh

#SGE_TASK_ID=22

# directories
#readonly in_dir_phased="data/phased"
#readonly in_dir_unphased="data/unphased/unfiltered"
#readonly vep_dir="data/vep/full"
#readonly spark_dir="data/tmp/spark"
readonly out_dir="data/saige/grm/input"

# hail script
#readonly hail_script="utils/hail_plink_export.py"

# input path
#readonly chr=${SGE_TASK_ID}
#readonly in_phased="${in_dir_phased}/ukb_wes_200k_phased_chr${chr}.1of1.vcf.gz"
#readonly in_unphased="${in_dir_unphased}/ukb_wes_200k_filtered_chr${chr}.mt"
#readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"

# output path
# path ukb_imp_eur_chr1_22_sparse_markers"
readonly out_prefix="${out_dir}/ukb_imp_eur_chr1_22_sparse_markers"
#readonly out="${out_prefix}.plink"

# SAIGE step paths
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"
#readonly step1_fitNULLGLMM="/well/lindgren/flassen/software/dev/SAIGE/step1_fitNULLGLMM.R"
#readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

# setup hail

# SECONDS=0
# set_up_hail
# mkdir -p ${out_dir}
# python3 "${hail_script}" \
#    --chrom ${chr} \
#    --input_path ${in_unphased} \
#    --input_type "mt" \
#    --get_europeans \
#    --missing 0.05 \
#    --out_prefix ${out_prefix} \
#    --out_type "plink"


set_up_RSAIGE
print_update "Generating GRM from plink files.. "
Rscript "${createSparseGRM}" \
	--plinkFile=${out_prefix} \
	--nThreads=19 \
	--outputPrefix=${out_prefix}	\
	--numRandomMarkerforSparseKin=1000	\
	--relatednessCutoff=0.125


print_update "Finished running HAIL for chr${chr}" "${SECONDS}"

