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
readonly grm_mtx=${8?Error: Missing arg8 (grm_mtx)}
readonly grm_sam=${9?Error: Missing arg9 (grm_sam)}
readonly phenotype=${10?Error: Missing arg9 (grm_sam)}

readonly csi="${vcf}.csi"
readonly step2_SPAtests="utils/saige/step2_SPAtests.R"
readonly shell_spa="scripts/conditional/common/_chr_spa.sh"
readonly rscript="scripts/conditional/common/03_spa_conditional.R"
readonly helper="scripts/conditional/common/_get_marker_and_pval.R"

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
  rm -f ${out}
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
   for chr in ${CHROMS}; do
      local out_spa="${spa_prefix}${chr}"
      local spa_copy="${spa_prefix}${chr}.txt"
      if [ ! -f "${out_spa}" ]; then
        Rscript "${step2_SPAtests}"  \
           --vcfFile="${vcf}" \
           --vcfFileIndex="${csi}" \
           --vcfField="GT" \
           --sparseGRMFile="${grm_mtx}" \
           --sparseGRMSampleIDFile="${grm_sam}"  \
           --chrom="chr${chr}" \
           --minMAF=0.0000001 \
           --minMAC="${min_mac}" \
           --GMMATmodelFile="${in_gmat}" \
           --varianceRatioFile="${in_var}" \
           --SAIGEOutputFile="${out_spa}" \
           --LOCO=FALSE \
           ${markers:+--condition "$markers"}
        # create backup of the chromosome, the other
        # is used for merging and then deleted subsequently
      else
        >&2 echo "${out_spa} already exists. Skipping.."
      fi
      cp "${out_spa}" "${spa_copy}"
   done
}



conditional_analysis() {
  
  # setup variables
  local i=0 
  local marker_list=""
  local current_marker=""
  local current_p=""

  >&2 echo "Beginning conditional analysis for ${phenotype}.."

  while [ $i -lt $max_iter ]; do

      true $(( i++ ))
      >&2 echo "[Iteration ${i}]: Current list of condtioning marker(s): ${marker_list}"

      local out_prefix_iter="${out_prefix}_i${i}"
      local out_mrg="${out_prefix_iter}.mrg"
      local out_r="${out_prefix_iter}.txt"

      # run SAIGE and merge output
      spa_chr_loop "${marker_list}" "${out_prefix_iter}"
      merge_spa_by_chr "${out_prefix_iter}_chrCHR.txt" "${out_mrg}"

      # set old values for comparison
      local old_marker=${current_marker}
      local old_p=${current_p}
    
      # format file so that the right P-value is always extracted 
      # regardless of conditional analysis being enabled or not 
      Rscript "${helper}" --spa_file "${out_mrg}" --p_cutoff "${P_cutoff}"  --out_file "${out_r}"
      local current_p=$( tail -n1 "${out_r}" | cut -f2 )
      local current_marker=$( tail -n1 "${out_r}" | cut -f1 )

      if [ ! -z "${current_marker}" ]; then
        if [ "${current_marker}" != "${old_marker}" ]; then
          if [ -z "${marker_list}" ]; then
            local marker_list="${current_marker}"
          else
            local marker_list="${marker_list},${current_marker}"
          fi
          >&2 echo "[i${i}]: New marker found '${current_marker} (P-value=${current_p}). Conditioning on ${marker_list} for ${phenotype}'"
        else
          >&2 echo "[i${i}]: No other markers passing sig threshold (P-value<${P_cutoff}). Ended loop with markers: ${marker_list} for ${phenotype}."
          break
        fi
      else
       >&2 echo "[i${i}]: No markers found at P-value<${P_cutoff}. Loop ended with markers ${marker_list} for ${phenotype}."
       break
      fi
      # clean up
      #rm -f "${out_prefix_iter}_chr"*
      #rm -f ${out_mrg}
    echo -e "${i}\t${current_marker}\t${current_p}\t${P_cutoff}" >> ${final_markers}
  done
}

###############
# main script #
##############

# setup permenant variables
readonly CHROMS=$(extract_chr_from_vcf ${vcf}) 
readonly final_markers="${out_prefix}.markers"
rm -f ${final_markers}

# setup saige
module purge
set_up_RSAIGE
set +eu

# Run conditional analysis
conditional_analysis ${vcf} ${out_prefix}


