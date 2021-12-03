#!/usr/bin/env bash
#
#
#$ -N spa_conditional
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spa_conditional.log
#$ -e logs/spa_conditional.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
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

# saige null model
readonly in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
readonly in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"

# in/out files 
#readonly vcf="${in_dir}/211111_genetic_intervals_${phenotype}.vcf.bgz"
#readonly out_prefix="${out_dir}/211111_significant_variants_${phenotype}"
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

# loop params

spa_loop() {
   printf 'VCF %s' "$1"
   printf 'Markers %s' "$2"
   printf 'Prefix %s' "$3"
   for chr in {1..22}; do
      local prefix_chr="${3}_chr${chr}"
      set -x
      Rscript "${step2_SPAtests}"  \
        --vcfFile=${1} \
        --vcfFileIndex="${1}.csi" \
        --vcfField="GT" \
        --chrom="chr${chr}" \
        --minMAF=0.0000001 \
        --minMAC=1 \
        --GMMATmodelFile=${in_gmat} \
        --varianceRatioFile=${in_var} \
        --SAIGEOutputFile=${3} \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE \
        --IsOutputNinCaseCtrl=TRUE \
        --IsOutputHetHomCountsinCaseCtrl=TRUE \
        --LOCO=FALSE \
        ${2:+--condition "$2"}
      set +x
   done 
}

spa_conditional() {
  vcf=${1}
  out_prefix=${2}
  local i=0
  local max=5
  local marker_list=""
  while [ $i -lt $max ]; do
      
      printf 'Iteration: %s' "$i"
      printf 'Markers: %s' "$marker_list"
      true $(( i++ ))
      
      local prefix_iter="${out_prefix}_i${i}"
      local prefix_iter_chr="prefix_iter_chr"
      local markers_combined="${prefix_iter}.txt"

      spa_loop ${vcf} ${marker_list} ${prefix_iter}
      cat ${prefix_iter_chr}* > ${markers_combined}
      rm  ${prefix_iter_chr}*
      
      marker=$( cat "${markers_combined}" | \
        awk -v P="${P_cutoff}" '$13 < P' | \
        sort -k 13,13 | \
        cut -d" " -f3 | \
        head -n1)

      if [ -z ${marker} ]; then
          break    
      fi    
      
      if [ -z ${marker_list} ]; then
          markers=${marker}
      else
          marker_list+=",${marker}"        
      fi        
      
  done
}

# loop params
readonly P_cutoff=0.0001
set_up_RSAIGE
mkdir -p ${out_dir}

annotation="synonymous"
vcf="${in_dir}/211111_intervals_${annotation}_${phenotype}.vcf.bgz"
out_prefix="${out_dir}/211111_spa_conditional_${annotation}_${phenotype}"
spa_conditional ${vcf} ${out_prefix}





