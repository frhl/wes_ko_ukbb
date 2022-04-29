#!/usr/bin/env bash
#
#
#$ -N _gene_spa
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_gene_spa.log
#$ -e logs/_gene_spa.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly spark_dir="data/tmp/spark"
readonly step2_SPAtests="utils/saige/step2_SPAtests.R"

readonly chr=${1?Error: Missing arg1 (chr)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly out_gene=${3?Error: Missing arg3 (out_gene)}
readonly out_spa_success=${4?Error: Missing arg3 (out_gene)}
readonly in_gmat=${5?Error: Missing arg4 (in_gmat)}
readonly in_var=${6?Error: Missing arg5 (in_var)}
readonly phenotype=${7?Error: Missing arg6 (phenotype)}
readonly gene=${8?Error: Missing arg7 (gene)}
readonly min_mac=${9?Error: Missing arg9 (min_mac)}

readonly var_bytes=$( file_size ${in_var} )
readonly gmat_bytes=$( file_size ${in_gmat} )

readonly id=${SGE_TASK_ID}
readonly vcf="${in_vcf}_${id}.vcf.gz"
readonly csi="${vcf}_${id}.csi"
readonly out_gene_task="${out_gene}_${id}.txt"
readonly out_file_success="${out_spa_success}_${id}.SUCCESS"
readonly out_file_failure="${out_spa_success}_${id}.FAILURE"

if [ ! -f ${out_gene_task} ]; then
  echo "var_bytes=${var_bytes} at ${in_var}"
  echo "gmat_bytes=${gmat_bytes} at ${in_gmat}"
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
    gzip ${out_gene_task}
  else
    touch ${out_file_failure}
  fi
else
  >&2 echo "${out_gene_task} already exists. Skipping.."
fi


