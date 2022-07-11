#!/usr/bin/env bash
#
#
#$ -N _spa_cond_rare
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_spa_cond_rare3.log
#$ -e logs/_spa_cond_rare3.errors.log
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
readonly grm_mtx=${6?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${7?Error: Missing arg7 (grm_sam)}
readonly min_mac=${8?Error: Missing arg8 (min_mac)} 
readonly in_markers_ac=${9?Error: Missing arg9 (markers_ac)} 
readonly out_prefix=${10?Error: Missing arg10 (out_prefix)}
readonly cond="${11}"
readonly cond_cat="${12}"
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
readonly markers_ac=$(echo ${in_markers_ac} | sed -e "s/CHR/${chr}/g")

readonly threads=$(( ${NSLOTS}-1 ))
readonly step2_SPAtests="utils/saige/step2_SPAtests_cond.R"
readonly rscript="scripts/conditional/rare/_variants_with_ac.R"

# condiitonal markers based on rare variants
readonly cond_chr=$(echo ${cond} | sed -e "s/CHR/${chr}/g")
>&2 echo ${cond_chr}
>&2 echo ${cond_cat}

# subset first by consequence
readonly markers_raw=$(zcat ${cond_chr} | grep -E "${cond_cat}" | cut -f3)
readonly markers_n=$(zcat ${cond_chr} | grep -E "${cond_cat}" | wc -l)
readonly markers_file="${out_prefix/CHR/${chr}}.markers"
echo ${markers_raw} > "${markers_file}"

# subset by non-monomorphic markers
set_up_rpy
readonly markers_pheno_file="${out_prefix/CHR/${chr}}.${phenotype}.markers"
Rscript "${rscript}" \
  --phenotype ${phenotype} \
  --current_marker_file ${markers_file} \
  --allele_count_file ${markers_ac} \
  --outfile ${markers_pheno_file} \
  --allele_count_threshold 2 #{min_mac}

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

rm ${markers_file}


