#!/usr/bin/env bash
#
#$ -N _fit_grm
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_fit_grm.log
#$ -e logs/_fit_grm.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 6
#$ -q short.qc@@short.hge
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly createSparseGRM="utils/saige/createSparseGRM.R"

readonly plink_file=${1?Error: Missing arg2 (in_vcf)}
readonly out_prefix=${3?Error: Missing arg2 (in_vcf)}
readonly threads="42"

mkdir -p ${out_dir}

SECONDS=0
set_up_RSAIGE
Rscript "${createSparseGRM}" \
    --plinkFile=${plink_file} \
    --nThreads=${threads} \
    --outputPrefix=${out_prefix} \
    --numRandomMarkerforSparseKin=2000 \
    --relatednessCutoff=0.05 \
    && print_update "Finished fitting GRM" ${SECONDS} \
    || raise_error "Fitting GRM failed"



