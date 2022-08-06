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

readonly threads=$(( ${NSLOTS}-1 ))
readonly step2_SPAtests="utils/saige/step2_SPAtests_cond.R"
readonly rscript_rare="scripts/conditional/rare/_spa_cond_rare.R"
readonly rscript_common="scripts/conditional/common/_spa_cond_common.R"
readonly rscript_merge="scripts/conditional/combined/_brute_force_cond.R"

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg3 (in_csi)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)}
readonly in_var=${5?Error: Missing arg5 (in_var)}
readonly grm_mtx=${6?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${7?Error: Missing arg7 (grm_sam)}
readonly min_mac=${8?Error: Missing arg8 (min_mac)}
readonly out_prefix=${9?Error: Missing arg9 (out_prefix)}
readonly in_markers_rare_ac=${10?Error: Missing arg10 (markers_rare_ac)}
readonly in_markers_rare_hash=${11?Error: Missing arg10 (markers_rare_ac)}
readonly markers_rare_by_gene=${12?Error: Missing arg10 (markers_rare_ac)}
readonly markers_rare_cond_min_mac=${13?Error: Missing arg10 (markers_rare_ac)}
readonly cond_rare_file=${14?Error: Missing arg11 (cond_rare_file)}
readonly cond_common_file=${15?Error: Missing arg12 (cond_common_file)}
readonly cond_annotation=${16?Error: Missing arg13 (cond_annotation)}
readonly chr=${SGE_TASK_ID}

# Need to change CHR input depending on current task-id
readonly gmat=$(echo ${in_gmat} | sed -e "s/CHR/${chr}/g")
readonly var=$(echo ${in_var} | sed -e "s/CHR/${chr}/g")
readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")
readonly markers_rare_ac=$(echo ${in_markers_rare_ac} | sed -e "s/CHR/${chr}/g")
readonly markers_rare_hash=$(echo ${in_markers_rare_hash} | sed -e "s/CHR/${chr}/g")

# Check that SAIGE step 1 has been run
readonly var_bytes=$( file_size ${var} )
readonly gmat_bytes=$( file_size ${gmat} )
>&2 echo "Found ${gmat} and ${var} with ${phenotype}"

# set up R
set_up_rpy

# create subset of rare (conding) markers to be used in analysis
readonly out_rare_markers_file="${out_prefix/CHR/${chr}}.rare.markers"
Rscript "${rscript_rare}" \
  --chromosome "chr${chr}" \
  --phenotype "${phenotype}" \
  --annotation "${cond_annotation}" \
  --path_ac_by_phenotypes "${markers_rare_ac}" \
  --path_hash_by_phenotypes "${markers_rare_hash}" \
  --path_markers_by_gene "${markers_rare_by_gene}" \
  --outfile "${out_rare_markers_file}" \
  --min_mac ${markers_rare_cond_min_mac}

# create a subset of common (non-coding) makers to be used in analysis
readonly out_common_markers_file="${out_prefix/CHR/${chr}}.common.markers"
readonly cond_chr_common_file=$(echo ${cond_common_file} | sed -e "s/CHR/${chr}/g")
Rscript ${rscript_common} \
  --infile ${cond_chr_common_file} \
  --phenotype ${phenotype} \
  --outfile ${out_common_markers_file} \
  --pheno_col 6

# merge the two subsets
readonly out_markers_file="${out_prefix/CHR/${chr}}.final.markers"
Rscript ${rscript_merge} \
  --file_rare_markers ${out_rare_markers_file} \
  --file_common_markers ${out_common_markers_file} \
  --outfile ${out_markers_file}


spa_test() {
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
       --condition_file "${out_markers_file}" \
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
  #spa_test
else
  >&2 echo "${out} already exists. Skipping.."
fi


