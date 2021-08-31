#!/usr/bin/env bash
#
#
#$ -N saige_hail
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o conda.log
#$ -e conda.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 22

#set -o errexit
#set -o nounset

source utils/bash_utils.sh

chr=22

readonly out_dir="data/saige/input"
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"
readonly out="${out_prefix}.plink"

# SAIGE step paths
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"

module load Anaconda3/2020.07
module load java/1.8.0_latest
source "/apps/eb/skylake/software/Anaconda3/2020.07/etc/profile.d/conda.sh"
conda activate RSAIGE


print_update "Generating GRM from plink files.. "
Rscript "${createSparseGRM}" \
	--plinkFile=${out_prefix} \
	--nThreads=1 \
	--outputPrefix=${out_prefix}	\
	--numRandomMarkerforSparseKin=1000	\
	--relatednessCutoff=0.125


print_update "Finished running HAIL for chr${chr}" "${SECONDS}"




