#!/usr/bin/env bash

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg3 (in_csi)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)} 
readonly in_var=${5?Error: Missing arg5 (in_var)} 
readonly grm_mtx=${6?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${7?Error: Missing arg7 (grm_sam)}
readonly min_mac=${8?Error: Missing arg8 (min_mac)} 
readonly in_group=${9?Error: Missing arg9 (path prefix for group file)}
readonly out_prefix=${10?Error: Missing arg10 (path prefix for saige output)}
readonly in_markers=${11} # optional conditioning markers
readonly chr=${SLURM_ARRAY_TASK_ID}

# chromosome specified at in_gmat and in_var
# when off chromosome PRS will beused
readonly gmat=$(echo ${in_gmat} | sed -e "s/CHR/${chr}/g")
readonly var=$(echo ${in_var} | sed -e "s/CHR/${chr}/g")

>&2 echo "${gmat} and ${var} with ${phenotype}"

readonly var_bytes=$( file_size ${var} )
readonly gmat_bytes=$( file_size ${gmat} )

readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly group=$(echo ${in_group} | sed -e "s/CHR/${chr}/g")
readonly markers=$(cat $(echo ${in_markers} | sed -e "s/CHR/${chr}/g"))
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

readonly threads=$(( ${NSLOTS}-1 ))
readonly step2_SPAtests="utils/saige/step2_SPAtests.R"

spa_set_test() {
   >&2 echo "var_bytes=${var_bytes} at ${var}"
   >&2 echo "gmat_bytes=${gmat_bytes} at ${gmat}"
   if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then 
     SECONDS=0
     set -x
     Rscript "${step2_SPAtests}"	\
       --vcfFile=${vcf} \
       --vcfFileIndex=${csi} \
       --vcfField="GT" \
       --sparseGRMFile=${grm_mtx} \
       --sparseGRMSampleIDFile=${grm_sam}  \
       --chrom="chr${chr}" \
       --minMAF=0 \
       --minMAC=${min_mac} \
       --GMMATmodelFile=${gmat} \
       --varianceRatioFile=${var} \
       --SAIGEOutputFile=${out} \
       --LOCO=FALSE \
       --groupFile=${group} \
       --annotation_in_groupTest=pLoF,damaging_missense,synonymous,damaging_missense:pLoF,synonymous:damaging_missense:pLoF \
       --maxMAF_in_groupTest=0.001,0.01,0.05
       ${markers:+--condition "$markers"} 
     set +x
   else
     raise_error "${var} or ${gmat} does not contain any bytes!"
   fi
}

if [ ! -f ${out} ]; then
   #set_up_RSAIGE 1.0.0
   set_up_RSAIGE
   spa_set_test
else
  >&2 echo "${out} already exists. Skipping.."
fi 


