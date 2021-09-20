#!/usr/bin/env bash
#
#
#$ -N mk_sparse_grm
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/create_sparse_grm.log
#$ -e logs/create_sparse_grm.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q short.qe

set -o errexit
set -o nounset

module purge
source utils/qsub_utils.sh
source utils/hail_utils.sh

# paths for hail script
readonly out_dir="data/saige/grm/input"
readonly out_prefix="${out_dir}/ukb_imp_eur_chr1_22_sparse_markers"
readonly hail_script="utils/create_grm_input.py"
readonly spark_dir="data/tmp/spark"
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"

# Generate a sequence of chromosomes to be included
chroms=$( seq 1 22 | tr '\n' ' ' )

# combine markers from UKBB imputed data
set_up_hail
set_up_pythonpath
mkdir -p ${out_dir}
python "${hail_script}" \
	--chroms ${chroms}  \
	--out_prefix ${out_prefix} \
	--subset_markers_by_kinship \
    --subset_samples_by_eur

print_update "Successfully combined .bgen files for GRM input." "${SECONDS}"

# generate GRM using SAIGE
set_up_RSAIGE
print_update "Generating GRM from plink files.. "
Rscript "${createSparseGRM}" \
	--plinkFile=${out_prefix} \
	--nThreads=4 \
	--outputPrefix=${out_prefix}	\
	--numRandomMarkerforSparseKin=1000	\
	--relatednessCutoff=0.125

