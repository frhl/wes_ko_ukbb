#!/usr/bin/env bash
#
#
#$ -N _array_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_array_permute.log
#$ -e logs/_array_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly chr=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg1 (phenotype)}
readonly in_tsv=${3?Error: Missing arg2 (in_vcf)}
readonly in_spa=${4?Error: Missing arg3 (in_csi)}
readonly out_prefix=${5?Error: Missing arg6 (path prefix for saige output)}

readonly threads=$(( ${NSLOTS}-1 ))
readonly step2_SPAtests="utils/saige/step2_SPAtests.R"


calc_permutations_required() 
{
  set_up_rpy
  Rscript ${rscript} \
    --in_spa ${in_spa} \
    --in_tsv ${in_tsv} \
    --out_prefix ${out_prefix} \
}




submit_permute_chunk() {

}






