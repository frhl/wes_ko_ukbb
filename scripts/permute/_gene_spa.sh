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

readonly chr=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg1 (phenotype)}
readonly out_gene=${3?Error: Missing arg1 (phenotype)}
readonly in_gmat=${4?Error: Missing arg2 (in_vcf)}
readonly in_var=${5?Error: Missing arg6 (path prefix for saige output)}
readonly phenotype=${6?Error: Missing arg6 (path prefix for saige output)}
readonly gene=${7?Error: Missing arg6 (path prefix for saige output)}
readonly n_tasks=${8?Error: Missing arg6 (path prefix for saige output)}
readonly min_mac=${9?Error: Missing arg6 (path prefix for saige output)}

readonly id=${SGE_TASK_ID}
readonly vcf="${in_vcf}_${id}of${n_tasks}.vcf.gz"
readonly csi="${vcf}.csi"
readonly out_gene_task="${out_gene}_${id}of${n_tasks}.txt"

set +eu
set_up_RSAIGE
set -eu

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

rm "${out_gene_task}.index"



