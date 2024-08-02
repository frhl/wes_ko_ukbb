#!/usr/bin/env bash
#

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly cluster=$( get_current_cluster)
readonly index=$( get_array_task_id )
readonly chr=$( get_chr ${index} )
readonly threads=( get_threads )

readonly step2_SPAtests="utils/saige/step2_SPAtests_cond.R"
readonly rscript="scripts/conditional/rare/_spa_cond_rare.R"

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg3 (in_csi)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)} 
readonly in_var=${5?Error: Missing arg5 (in_var)} 
readonly grm_mtx=${6?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${7?Error: Missing arg7 (grm_sam)}
readonly min_mac=${8?Error: Missing arg8 (min_mac)} 
readonly in_markers_ac=${9?Error: Missing arg9 (markers_ac)} 
readonly in_markers_hash=${10?Error: Missing arg9 (markers_hash)} 
readonly markers_by_gene=${11?Error: Missing arg9 (markers_hash)} 
readonly markers_cond_min_mac=${12?Error: Missing arg9 (markers_hash)} 
readonly out_prefix=${13?Error: Missing arg10 (out_prefix)}
readonly cond_markers="${14}"
readonly cond_annotation="${15}"



# Change CHR to current chromosome based on task-ID
readonly gmat=$(echo ${in_gmat} | sed -e "s/CHR/${chr}/g")
readonly var=$(echo ${in_var} | sed -e "s/CHR/${chr}/g")
readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")
readonly markers_ac=$(echo ${in_markers_ac} | sed -e "s/CHR/${chr}/g")
readonly markers_hash=$(echo ${in_markers_hash} | sed -e "s/CHR/${chr}/g")

# check that saige step 1 exists
readonly var_bytes=$( file_size ${var} )
readonly gmat_bytes=$( file_size ${gmat} )

# setup R
set_up_rpy

# todo: Only start conditional analysis if markers_by_gene has been generated (which
# is currently only generated for markers that are significant in main analysis)

# subset variants to be used for conditional analysis
readonly cond_markers_chr=$(echo ${cond_markers} | sed -e "s/CHR/${chr}/g")
readonly markers_pheno_file="${out_prefix/CHR/${chr}}.rare.markers"
Rscript "${rscript}" \
  --chromosome "chr${chr}" \
  --phenotype "${phenotype}" \
  --annotation "${cond_annotation}" \
  --path_markers_in_chrom "${cond_markers_chr}" \
  --path_ac_by_phenotypes "${markers_ac}" \
  --path_hash_by_phenotypes "${markers_hash}" \
  --path_markers_by_gene "${markers_by_gene}" \
  --outfile "${markers_pheno_file}" \
  --min_mac "${markers_cond_min_mac}"

spa_test() {
  echo "var_bytes=${var_bytes} at ${var}"
  echo "gmat_bytes=${gmat_bytes} at ${gmat}"
  if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then 
    SECONDS=0
    Rscript "${step2_SPAtests}" \
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
       --LOCO=FALSE \
       --condition_file "${markers_pheno_file}" \
       && print_update "Finished saddle-point approximation for chr${chr}" ${SECONDS} \
       || raise_error "Saddle-point approximation for chr${chr} failed"
  else
    raise_error "${var} or ${gmat} does not contain any bytes!"
  fi
}
if [ ! -f ${out} ]; then
  set +eu
  conda deactivate
  set_up_RSAIGE
  set -eu
  spa_test
else
  >&2 echo "${out} already exists. Skipping.."
fi 



