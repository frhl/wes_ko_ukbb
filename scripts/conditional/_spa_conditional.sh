#!/usr/bin/env bash
#
#
#$ -N _spa_conditional
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_conditional.log
#$ -e logs/spa_conditional.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

source utils/bash_utils.sh

#readonly in_dir=${1?Error: Missing arg1 (in dir)}
#readonly out_dir=${2?Error: Missing arg2 (out_dir)}
#readonly pheno_dir=${3?Error: Missing arg3 (pheno_dir)}
#readonly step1_dir=${4?Error: Missing arg4 (step1_dir)}
#readonly pheno_list=${1?Error: Missing arg1 (pheno_dir)}
readonly in_gmat=${1?Error: Missing arg2 (in_gmat)}
readonly in_var=${2?Error: Missing arg3 (in_var)}
readonly vcf=${3?Error: Missing arg4 (vcf)}
readonly out_prefix=${4?Error: Missing arg5 (out_prefix)}
readonly P_cutoff=${5?Error: Missing arg6 (P_cutoff)}
readonly max_iter=${6?Error: Missing arg7 (max_iter)}

readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"
readonly rscript="scripts/conditional/03_spa_conditional.R"

spa_chr_loop() {
   for chr in {1..22}; do
      local prefix_chr="${3}_chr${chr}"
      set -x
      Rscript "${step2_SPAtests}"  \
        --vcfFile=${1} \
        --vcfFileIndex="${1}.csi" \
        --vcfField="GT" \
        --chrom="chr${chr}" \
        --minMAF=0.000001 \
        --minMAC=1 \
        --GMMATmodelFile=${in_gmat} \
        --varianceRatioFile=${in_var} \
        --SAIGEOutputFile=${prefix_chr} \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE \
        --IsOutputNinCaseCtrl=TRUE \
        --IsOutputHetHomCountsinCaseCtrl=TRUE \
        --LOCO=FALSE \
        ${2:+--condition "$2"}
      set +x
   done
   >&2 echo "Markers: ${2}"
   >&2 echo "VCF: ${1}"
   >&2 echo "Prefix: ${3}"
}

conditional_analysis() {
  markers_conditional="${out_prefix}_conditioning_markers.txt"
  rm ${markers_conditional}

  i=0 #
  marker_list="" # markers for SAIGE
  while [ $i -lt $max_iter ]; do

      true $(( i++ ))
      >&2 echo "[Iteration ${i}]: Beginning new iteration.."
      >&2 echo "[Iteration ${i}]: Current condtioning marker: ${current_marker}"
      >&2 echo "[Iteration ${i}]: Current condtioning marker list: ${marker_list}"

      # setup prefixes for temporary files
      prefix_iter="${out_prefix}_i${i}"
      prefix_iter_chr="${prefix_iter}_chr"
      markers_combined="${prefix_iter}.txt"

      # Run saddle point approximation using condtioning markers
      spa_chr_loop "${vcf}" "${marker_list}" "${prefix_iter}"
      Rscript "${rscript}" --prefix ${prefix_iter}
      rm  "${prefix_iter_chr}"*

      # Save top marker within P-value threshold
      cat "${markers_combined}" | \
        awk -v P="${P_cutoff}" '$22 < P' | \
        head -n1 >> "${markers_conditional}"

      # select current marker
      old_marker=${current_marker}
      current_marker=$( tail -n1 ${markers_conditional} | awk '{print $1":"$2"_"$4"/"$5}')
      current_p=$( tail -n1 ${markers_conditional} | cut -d" " -f22)

      if [[ "${current_marker}" != "${old_marker}" ]]; then
        if [[ "${marker_list}" == "" ]]; then
          marker_list="${current_marker}"
        else
          marker_list="${marker_list},${current_marker}"
        fi
        >&2 echo "[Iteration ${i}]: New marker found '${current_marker} (P-value=${current_p}). Conditioning on ${marker_list}'"
      else
        >&2 echo "[Iteration ${i}]: No other markers passing sig threshold (P-value<${P_cutoff}). Ended loop with markers: ${marker_list}"
        break
      fi

  done
}

# Run analysis
set_up_RSAIGE
conditional_analysis ${vcf} ${out_prefix}

