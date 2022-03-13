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

readonly threads=$(( ${NSLOTS}-1 ))
readonly step2_SPAtests="utils/saige/step2_SPAtests.R"
readonly rscript="scripts/conditional/03_spa_conditional.R"

spa_chr_loop() {
   local cond_vcf=${1}
   local cond_markers=${2}
   local cond_csi="${vcf}.csi"
   for chr in {1..22}; do
      SECONDS=0
      local prefix_chr="${3}_chr${chr}"
      Rscript "${step2_SPAtests}"  \
         --vcfFile=${cond_vcf} \
         --vcfFileIndex=${cond_csi} \
         --vcfField="GT" \
         --chrom="chr${chr}" \
         --minMAF=0.0000001 \
         --minMAC=${min_mac} \
         --GMMATmodelFile=${in_gmat} \
         --varianceRatioFile=${in_var} \
         --SAIGEOutputFile=${prefix_chr} \
         --LOCO=FALSE \
         ${cond_markers:+--condition "$cond_markers"} \
         && print_update "Finished saddle-point approximation for chr${chr}" ${SECONDS} \
         || raise_error "Saddle-point approximation for chr${chr} failed" 
   done
   >&2 echo "Markers: ${2}"
   >&2 echo "VCF: ${1}"
   >&2 echo "Prefix: ${3}"
}

conditional_analysis() {
  local markers_conditional="${out_prefix}_conditioning_markers.txt"
  rm -f ${markers_conditional}

  local i=0 
  local marker_list=""
  local current_marker=""
  while [ $i -lt $max_iter ]; do

      true $(( i++ ))
      >&2 echo "[Iteration ${i}]: Beginning new iteration.."
      >&2 echo "[Iteration ${i}]: Current condtioning marker: ${current_marker}"
      >&2 echo "[Iteration ${i}]: Current condtioning marker list: ${marker_list}"

      # setup prefixes for temporary files
      local prefix_iter="${out_prefix}_i${i}"
      local prefix_iter_chr="${prefix_iter}_chr"
      local markers_combined="${prefix_iter}.txt"

      # Run saddle point approximation using condtioning markers
      spa_chr_loop "${vcf}" "${marker_list}" "${prefix_iter}"
      Rscript "${rscript}" --prefix ${prefix_iter}
      rm  "${prefix_iter_chr}"*

      # Save top marker within P-value threshold
      cat "${markers_combined}" | \
        awk -v P="${P_cutoff}" '$22 < P' | \
        head -n1 >> "${markers_conditional}"

      # select current marker
      local old_marker=${current_marker}
      local current_marker=$( tail -n1 ${markers_conditional} | awk '{print $1":"$2"_"$4"/"$5}')
      local current_p=$( tail -n1 ${markers_conditional} | cut -d" " -f22)

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

  done
}

# Run analysis
set +eu
set_up_RSAIGE
set -eu
conditional_analysis ${vcf} ${out_prefix}

