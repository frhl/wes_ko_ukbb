#!/usr/bin/env bash
#
#
#$ -N create_grm
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/create_grm.log
#$ -e logs/create_grm.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe

#set -o errexit
#set -o nounset

module purge
source utils/bash_utils.sh

# paths for hail script
readonly out_dir="data/saige/grm/input"
readonly out_prefix="${out_dor}/ukb_imp_sparse_markers"
readonly hail_script="utils/create_grm_input.py"
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"

# Generate a sequence of chromosomes to be included
chroms=$( seq 21 22 | tr '\n' ' ' )

mkdir -p ${out_dir}
python3 "${hail_script}" \
    --chroms ${chroms}  \
    --out_prefix ${out_prefix} \
    --subset_markers_by_kinship


# compute GRM
#conda deacticate
#set_up_RSAIGE
#print_update "Generating GRM from plink files.. "
#Rscript "${createSparseGRM}" \
#	--plinkFile=${out_prefix} \
#	--nThreads=4 \
#	--outputPrefix=${out_prefix}	\
#	--numRandomMarkerforSparseKin=1000	\
#	--relatednessCutoff=0.125


print_update "Finished running HAIL for chr${chr}" "${SECONDS}"

