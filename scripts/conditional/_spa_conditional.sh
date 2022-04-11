#!/usr/bin/env bash
#
#
#$ -N _spa_conditional
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_spa_conditional.log
#$ -e logs/_spa_conditional.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly in_gmat=${1?Error: Missing arg1 (in_gmat)}
readonly in_var=${2?Error: Missing arg2 (in_var)}
readonly vcf=${3?Error: Missing arg3 (vcf)}
readonly out_prefix=${4?Error: Missing arg4 (out_prefix)}
readonly P_cutoff=${5?Error: Missing arg5 (P_cutoff)}
readonly max_iter=${6?Error: Missing arg6 (max_iter)}
readonly min_mac=${7?Error: Missing arg7 (min_mac)}

readonly csi="${vcf}.csi"
readonly step2_SPAtests="utils/saige/step2_SPAtests.R"
readonly shell_spa="scripts/conditional/_chr_spa.sh"
readonly rscript="scripts/conditional/03_spa_conditional.R"

# A function to extract all the (unique) chromsomes
extract_chr_from_vcf() {
  module purge
  module load BCFtools/1.12-GCC-10.3.0 
  bcftools query -f '%CHROM\n' ${1} | sort | uniq | sed -e "s/chr//g" |tr "\n" " "
}

# Merge saddle point approximation by chromosome
# expects chromosome to the only difference in filename
merge_spa_by_chr(){
  local prefix=${1?Error: expected prefix (arg1)}
  local out=${2?Error: expected out file (arg2)}
  local first_iter="Y"
  for chr in {1..22}; do
     local file=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")
     if [ -f ${file} ]; then
       if [ ${first_iter} == "Y" ]; then
          cat "${file}" | head -n 1  >> "${out}"
          local first_iter="N"
       fi
       cat "${file}" | grep -v "MarkerID"  >> "${out}"
       rm -f "${file}"
       rm -f "${file}.index"
     fi
  done
}

# perform saddlepoint approximation loop across all 
# chromosomes that are present in the current vcf
spa_chr_loop() {
   local markers="${1}"
   local spa_prefix="${2}_chr"
   set +eu
   for chr in ${CHROMS}; do
      if [ ! -f "${spa_prefix}${chr}" ]; then
        Rscript "${step2_SPAtests}"  \
           --vcfFile="${vcf}" \
           --vcfFileIndex="${csi}" \
           --vcfField="GT" \
           --chrom="chr${chr}" \
           --minMAF=0.0000001 \
           --minMAC="${min_mac}" \
           --GMMATmodelFile="${in_gmat}" \
           --varianceRatioFile="${in_var}" \
           --SAIGEOutputFile="${spa_prefix}${chr}" \
           --LOCO=FALSE \
           ${markers:+--condition "$markers"}
      fi
      cp "${spa_prefix}${chr}" "${spa_prefix}${chr}.txt"
   done
   set -eu
}

qsub_spa_chr_loop() {
# deprecated for now ..
   local markers="${1}"
   local spa_prefix="${2}_chr"
   qsub -N "_cond_spa_${markers}" \
        -t 1-22 \
        -q "short.qf" \
        -pe shmem 1 \
        "${shell_spa}" \
        "${vcf}" \
        "${in_gmat}" \
        "${in_var}" \
        "${phenotype}" \
        "${prefix}" \
        "${markers}"
}



conditional_analysis() {
  
  local i=0 
  local marker_list=""
  local current_marker=""
  local current_p=""

  while [ $i -lt $max_iter ]; do

      true $(( i++ ))
      >&2 echo "[Iteration ${i}]: Current list of condtioning marker(s): ${marker_list}"

      local out_prefix_iter="${out_prefix}_i${i}"
      local out_prefix_mrg="${out_prefix_iter}.mrg"

      spa_chr_loop "${marker_list}" "${out_prefix_iter}"
      merge_spa_by_chr "${out_prefix_iter}_chrCHR.txt" "${out_prefix_mrg}"

      local old_marker=${current_marker}
      local old_p=${current_p}

      if [ ${i} -eq 1 ]; then

        cat ${out_prefix_mrg} | awk -v P="${P_cutoff}" '$13 < P' | head -n1 >> ${markers_conditional}
        local current_p=$( tail -n1 ${markers_conditional} | cut -f13)
        local current_marker=$( tail -n1 ${markers_conditional} | awk '{print $1":"$2":"$4":"$5}')

      else

        cat ${out_prefix_mrg} | awk -v P="${P_cutoff}" '$18 < P' | head -n1 >> ${markers_conditional} 
        local current_p=$( tail -n1 ${markers_conditional} | cut -f18)
        local current_marker=$( tail -n1 ${markers_conditional} | awk '{print $1":"$2":"$4":"$5}')

      fi

      if [[ "${current_marker}" != "${old_marker}" ]]; then
        if [[ "${marker_list}" == "" ]]; then
          local marker_list="${current_marker}"
        else
          local marker_list="${marker_list},${current_marker}"
        fi
        >&2 echo "[Iteration ${i}]: New marker found '${current_marker} (P-value=${current_p}). Conditioning on ${marker_list}'"
      else
        >&2 echo "[Iteration ${i}]: No other markers passing sig threshold (P-value<${P_cutoff}). Ended loop with markers: ${marker_list}"
        break
      fi
      # clean up
      rm -f "${out_prefix_iter}_chr"*
      rm -f ${out_prefix_mrg}
  done
  echo ${marker_list} > "${out_prefix}.markers"
}

###############
# main script #
##############

# setup constant files 
readonly CHROMS=$(extract_chr_from_vcf ${vcf}) 
readonly markers_conditional="${out_prefix}.spa"
rm -f ${markers_conditional}

# setup saige
module purge
set_up_RSAIGE

# Run conditional analysis
conditional_analysis ${vcf} ${out_prefix}


