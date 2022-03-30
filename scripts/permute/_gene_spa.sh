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

readonly id=${SGE_TASK_ID}
readonly vcf="${in_vcf}_${id}.vcf.gz"
readonly csi="${vcf}_${id}.csi"
readonly out_gene_task="${out_gene}_${id}.txt"
readonly out_file_success="${out_spa_success}_${id}.SUCCESS"

set_up_RSAIGE

if [ ! -f ${out_gene_task} ]; then
  SECONDS=0
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
  touch ${out_file_success}
else
  >&2 echo "${out_gene_task} already exists. Skipping.."
fi



