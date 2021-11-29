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


readonly out_dir=""
readonly out=""

# select phenotype (1-42)
readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

# get common variant files
readonly in_dir="data/conditional/common"
readonly vcf="${in_dir}/211111_genetic_intervals_${phenotype}.vcf.bgz"
readonly csi="${vcf}.csi"

# setup input phenotypes
readonly in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"

# conditional parameters
readonly max_iter=10
readonly P_cutoff=0.001

# SAIGE paths
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

spa_test() {
   SECONDS=0
   variants="${1}"
   iter="${2}"
   mv "${out_tmp}" "${out_tmp}_i${iter}"
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
     --condition "${variants}"
   set +x
   duration=${SECONDS}
   print_update "Finished SPA conditional on variants ${variants} at ${out_tmp} (${duration})"
}

get_spa_marker() {
  echo $( cat ${out_tmp} | \
    awk -v P="${P_cutoff}" '$13 < P' | \
    sort -k 13,13 | \
    cut -d" " -f3 | \
    head -n1)
}



iter=0
marker=$( get_spa_marker )
set_up_RSAIGE
while [ ! -z ${marker} ]; do  
  
  iter=`expr $iter + 1`
  
  # Keep track of variants
  if [ -z ${markers} ]; do
    markers="${marker}"
  else
    markers+=",${marker}"
  fi

  # perform SPA conditioning on current markers 
  spa_test "${markers}" "${iter}"

  # update marker
  marker=$( get_spa_marker )
 
  echo "${markers}" 
  if [ $iter -gt $max_iter ]; then
    echo "Exiting loop..."
    break
  fi
done 






