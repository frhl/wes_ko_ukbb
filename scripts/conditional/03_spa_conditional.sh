#!/usr/bin/env bash
#
#
#$ -N spa_conditional
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_conditional.log
#$ -e logs/spa_conditional.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge
#$ -t 3

source utils/bash_utils.sh

# Directories
readonly in_dir="data/conditional/common"
readonly out_dir="data/conditional/common"
readonly pheno_dir="data/phenotypes"
readonly step1_dir="data/saige/output/combined/binary/step1"

# setup phenotype
readonly pheno_list="${pheno_dir}/UKBB_WES200k_binary_phenotypes_header.txt"
readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

# saige parameters
readonly in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

spa_chr_loop() {
   for chr in {5..6}; do
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
  vcf=${1}
  out_prefix=${2}
  markers_conditional="${out_prefix}_conditioning_markers.txt"
  rm ${markers_conditional}
  
  i=0
  max=2
  marker_list=""
  while [ $i -lt $max ]; do
      
      true $(( i++ ))
      
      >&2 echo "Iteration: ${i}"
      >&2 echo "Condtioning marker: ${current_marker}"
      >&2 echo "Condtioning marker list: ${marker_list}"
      
      # setup prefixes for temporary files
      prefix_iter="${out_prefix}_i${i}"
      prefix_iter_chr="${prefix_iter}_chr"
      markers_combined="${prefix_iter}.txt"

      # Run saddle point approximation using condtioning markers
      spa_chr_loop "${vcf}" "${marker_list}" "${prefix_iter}"
      cat "${prefix_iter_chr}"* > ${markers_combined}
      rm  "${prefix_iter_chr}"*
      
      # Save top marker within P-value threshold
      cat "${markers_combined}" | \
        awk -v P="${P_cutoff}" '$13 < P' | \
        sort -k 13,13 | \
        head -n1 >> "${markers_conditional}" 

      # select current marker
      old_marker=${current_marker}
      current_marker=$( tail -n1 ${markers_conditional} | awk '{print $1":"$2"_"$4"/"$5}')
      >&2 echo "old_marker=${old_marker}"
      >&2 echo "current_marker=${current_marker}"
      if [[ "${current_marker}" != "${old_marker}" ]]; then
        marker_list=$( cat ${markers_conditional} | cut -d" " -f3 | tr "\n" "," | sed 's/,$//' )
        >&2 echo "New marker found '${current_marker}. Repeating loop with ${marker_list}'"
      else
        >&2 echo "Ended loop after ${i} iterations with markers: ${marker_list}"
        break    
      fi

  done
}

# loop params
set_up_RSAIGE
readonly P_cutoff=0.0005
mkdir -p ${out_dir}

annotation="synonymous"
vcf="${in_dir}/211111_intervals_${annotation}_${phenotype}.vcf.bgz"
out_prefix="${out_dir}/211111_spa_conditional_${annotation}_${phenotype}"
conditional_analysis ${vcf} ${out_prefix}





