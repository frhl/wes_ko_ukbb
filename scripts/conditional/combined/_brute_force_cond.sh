#!/usr/bin/env bash
#
#
#$ -N _brute_force_cond
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_brute_force_cond.log
#$ -e logs/_brute_force_cond.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg3 (in_csi)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)} 
readonly in_var=${5?Error: Missing arg5 (in_var)} 
readonly min_mac=${6?Error: Missing arg6 (min_mac)} 
readonly out_prefix=${7?Error: Missing arg7 (path prefix for saige output)}
readonly sorted_markers=${8?Error: Missing arg7 (path prefix for saige output)}
readonly cond="${9}"
readonly cond_cat="${10}"
readonly chr=${SGE_TASK_ID}


# chromosome specified at in_gmat and in_var
# when off chromosome PRS will beused
readonly gmat=$(echo ${in_gmat} | sed -e "s/CHR/${chr}/g")
readonly var=$(echo ${in_var} | sed -e "s/CHR/${chr}/g")

>&2 echo "${gmat} and ${var} with ${phenotype}"

readonly var_bytes=$( file_size ${var} )
readonly gmat_bytes=$( file_size ${gmat} )

readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

readonly threads=$(( ${NSLOTS}-1 ))
readonly step2_SPAtests="utils/saige/step2_SPAtests_cond.R"

# condiitonal markers based on rare variants
readonly cond_chr=$(echo ${cond} | sed -e "s/CHR/${chr}/g")
readonly markers_raw=$(zcat ${cond_chr} | grep -E "${cond_cat}" | cut -f3)
readonly markers_n=$(zcat ${cond_chr} | grep -E "${cond_cat}" | wc -l)
readonly markers_file="${out_prefix}.markers"


spa_test() {
  echo "var_bytes=${var_bytes} at ${var}"
  echo "gmat_bytes=${gmat_bytes} at ${gmat}"
  if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then 
    SECONDS=0
    Rscript "${step2_SPAtests}"	\
       --vcfFile=${vcf} \
       --vcfFileIndex=${csi} \
       --vcfField="DS" \
       --chrom="chr${chr}" \
       --minMAF=0.0000001 \
       --minMAC=${min_mac} \
       --GMMATmodelFile=${gmat} \
       --varianceRatioFile=${var} \
       --SAIGEOutputFile=${out} \
       --LOCO=FALSE \
       --condition_file "${markers_file}" \
       && print_update "Finished saddle-point approximation for chr${chr}" ${SECONDS} \
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

rm ${markers_file}


