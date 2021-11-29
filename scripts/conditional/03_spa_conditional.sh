#!/usr/bin/env bash
#
#
#$ -N spa_conditional
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_conditional.log
#$ -e logs/spa_conditional.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

source utils/bash_utils.sh
source utils/hail_utils.sh

# arguments
readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg2 (in_csi)}
readonly in_gmat=${4?Error: Missing arg3 (in_gmat)} 
readonly in_var=${5?Error: Missing arg4 (in_var)} 
readonly out=${6?Error: Missing arg5 (path prefix for saige output)}

readonly out_tmp="${out}_"

iter=0
max_iter=10
P_cutoff=0.001

# subsitute each chromosome
#readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
#readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")

# SAIGE paths
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

spa_test() {
   SECONDS=0
   print_update "SAIGE-SPA: Writing '${out}' using '${vcf}'"
   rm ${out_tmp}
   set -x
   Rscript "${step2_SPAtests}"  \
     --vcfFile=${vcf} \
     --vcfFileIndex=${csi} \
     --vcfField="DS" \
     --chrom="chr${chr}" \
     --minMAF=0.0000001 \
     --minMAC=3 \
     --GMMATmodelFile=${in_gmat} \
     --varianceRatioFile=${in_var} \
     --SAIGEOutputFile=${out_tmp} \
     --numLinesOutput=2 \
     --IsOutputAFinCaseCtrl=TRUE \
     --IsOutputNinCaseCtrl=TRUE \
     --IsOutputHetHomCountsinCaseCtrl=TRUE \
     --LOCO=FALSE \
     --condition "${1}"
   set +x
   duration=${SECONDS}
   print_update "Finished SPA test for variants ${1} with ${out} at ${duration}"
}

get_spa_marker() {
  echo $( cat ${out_tmp} | \
    awk -v P="${P_cutoff}" '$13 < P' | \
    sort -k 13,13 | \
    cut -d" " -f3 | \
    head -n1)
}


set_up_RSAIGE

markers=""
marker=$( get_spa_marker )
while [ ! -z ${marker} ]; do  
  markers+=",${marker}"
  spa_test "${markers}"
  marker=$( get_spa_marker )
  echo "${markers}" 
  iter=`expr $iter + 1`
  if [ $iter -gt $max_iter ]; then
    echo "Exiting loop..."
    break
  fi

done 






