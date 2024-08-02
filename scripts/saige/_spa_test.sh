#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg3 (in_csi)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)} 
readonly in_var=${5?Error: Missing arg5 (in_var)} 
readonly grm_mtx=${6?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${7?Error: Missing arg7 (grm_sam)}
readonly min_mac=${8?Error: Missing arg6 (min_mac)} 
readonly out_prefix=${9?Error: Missing arg7 (path prefix for saige output)}
readonly in_markers="${10}" # optional conditioning markers
readonly chr=$( get_array_task_id )


# chromosome specified at in_gmat and in_var
# when off chromosome PRS will beused
readonly gmat=$(echo ${in_gmat} | sed -e "s/CHR/${chr}/g")
readonly var=$(echo ${in_var} | sed -e "s/CHR/${chr}/g")

>&2 echo "${gmat} and ${var} with ${phenotype}"

readonly var_bytes=$( file_size ${var} )
readonly gmat_bytes=$( file_size ${gmat} )

readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly markers=$(cat $(echo ${in_markers} | sed -e "s/CHR/${chr}/g"))
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

readonly step2_SPAtests="utils/saige/step2_SPAtests.R"

spa_test() {
  echo "var_bytes=${var_bytes} at ${var}"
  echo "gmat_bytes=${gmat_bytes} at ${gmat}"
  if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then 
  SECONDS=0
  Rscript "${step2_SPAtests}"	\
     --vcfFile=${vcf} \
     --vcfFileIndex=${csi} \
     --vcfField="DS" \
     --sparseGRMFile=${grm_mtx} \
     --sparseGRMSampleIDFile=${grm_sam}  \
     --chrom="chr${chr}" \
     --minMAF=0.0000001 \
     --minMAC=${min_mac} \
     --GMMATmodelFile=${gmat} \
     --varianceRatioFile=${var} \
     --SAIGEOutputFile=${out} \
     --is_output_moreDetails=TRUE \
     --LOCO=FALSE\
     ${markers:+--condition "$markers"} \
    && print_update "Finished saddle-point approximation. Writing to ${out}" ${SECONDS} \
     || raise_error "Saddle-point approximation for chr${chr} failed"
  else
    raise_error "${var} or ${gmat} does not contain any bytes!"
  fi
}
 
if [ ! -f ${out} ]; then
   set_up_RSAIGE
   spa_test
else
  >&2 echo "${out} already exists. Skipping.."
fi 


