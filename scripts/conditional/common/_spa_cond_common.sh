#!/usr/bin/env bash
#
# Note: we assume that all phenotypes have common markers present, i.e.
# there is no need to subset to non-monomorphic markers in a phenotype
# dependent manner (see scripts/conditional/rare/_spa_cond_rare.sh).
#
#$ -N _spa_cond_rare
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_spa_cond_common.log
#$ -e logs/_spa_cond_common.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly threads=$(( ${NSLOTS}-1 ))
readonly step2_SPAtests="utils/saige/step2_SPAtests_cond.R"
readonly read_markers="scripts/conditional/utils/_read_markers.R"

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg3 (in_csi)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)} 
readonly in_var=${5?Error: Missing arg5 (in_var)} 
readonly grm_mtx=${6?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${7?Error: Missing arg7 (grm_sam)}
readonly min_mac=${8?Error: Missing arg8 (min_mac)} 
readonly out_prefix=${9?Error: Missing arg10 (out_prefix)}
readonly cond_markers="${10}"
readonly cond_cat="${11}"
readonly chr=${SGE_TASK_ID}

readonly gmat=$(echo ${in_gmat} | sed -e "s/CHR/${chr}/g")
readonly var=$(echo ${in_var} | sed -e "s/CHR/${chr}/g")

readonly var_bytes=$( file_size ${var} )
readonly gmat_bytes=$( file_size ${gmat} )

readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

readonly cond_markers_chr=$(echo ${cond_markers} | sed -e "s/CHR/${chr}/g") 

spa_test() {

  # generate a list of conditioning markers if available
  local out_markers="${out_prefix/CHR/${chr}}.markers"
  Rscript ${read_markers} \
    --infile ${cond_markers_chr} \
    --phenotype ${phenotype} \
    --outfile ${out_markers} \
    --pheno_col 5

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
       --condition_file ${out_markers} \
       && print_update "Finished saddle-point approximation for chr${chr}" ${SECONDS} \
       || raise_error "Saddle-point approximation for chr${chr} failed"
  else
    raise_error "${var} or ${gmat} does not contain any bytes!"
  fi
}
if [ ! -f ${out} ]; then
  set +eu
  set_up_RSAIGE
  set -eu
  spa_test
else
  >&2 echo "${out} already exists. Skipping.."
fi 

rm ${markers_file}


