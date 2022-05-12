#!/usr/bin/env bash
#
#
#$ -N _chr_spa
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_chr_spa.log
#$ -e logs/_chr_spa.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly step2_SPAtests="utils/saige/step2_SPAtests.R"

readonly in_vcf=${1?Error: Missing arg2 (in_vcf)}
readonly in_gmat=${2?Error: Missing arg4 (in_gmat)}
readonly in_var=${3?Error: Missing arg5 (in_var)}
readonly phenotype=${4?Error: Missing arg6 (phenotype)}
readonly min_mac=${5?Error: Missing arg9 (min_mac)}
readonly out_prefix=${6?Error: Missing arg9 (min_mac)}
readonly markers=${7}

readonly var_bytes=$( file_size ${in_var} )
readonly gmat_bytes=$( file_size ${in_gmat} )

readonly chr=${SGE_TASK_ID}
readonly vcf="${in_vcf}_${id}.vcf.gz"
readonly csi="${vcf}_${id}.csi"
readonly out_gene_task="${out_gene}_${id}.txt"
readonly out_file_success="${out_spa_success}_${id}.SUCCESS"
readonly out_file_failure="${out_spa_success}_${id}.FAILURE"

if [ ! -f ${out_gene_task} ]; then
  if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then
    SECONDS=0
    set_up_RSAIGE
    Rscript "${step2_SPAtests}"  \
       --vcfFile=${vcf} \
       --vcfFileIndex=${csi} \
       --vcfField="DS" \
       --chrom="chr${chr}" \
       --minMAF=0.000000001 \
       --minMAC=${min_mac} \
       --GMMATmodelFile=${in_gmat} \
       --varianceRatioFile=${in_var} \
       --SAIGEOutputFile=${out_gene_task} \
       --LOCO=FALSE\
       && print_update "Finished saddle-point approximation for chr${chr}" ${SECONDS} \
       || raise_error "Saddle-point approximation for chr${chr} failed"
    rm -f "${out_gene_task}.index"
  else
    touch ${out_file_failure}
  fi
else
  >&2 echo "${out_gene_task} already exists. Skipping.."
fi

