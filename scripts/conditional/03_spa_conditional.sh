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
readonly vcf="${in_dir}/211111_genetic_intervals_${phenotype}.vcf.bgz"
readonly csi="${vcf}.csi"
readonly out="${out_dir}/211111_significant_variants_${phenotype}.txt"
readonly out_prefix="${out_dir}/211111_significant_variants_${phenotype}"
readonly out_tmp="${out}.tmp"

# loop params
readonly P_cutoff=0.0001
readonly step2_SPAtests="/well/lindgren/flassen/software/dev/SAIGE/extdata/step2_SPAtests.R"

#spa_test() {
#  set -x
#  Rscript "${step2_SPAtests}"  \
#    --vcfFile=${vcf} \
#    --vcfFileIndex=${csi} \
#    --vcfField="GT" \
#    --chrom="chr${chr}" \
#    --minMAF=0.0000001 \
#    --minMAC=1 \
#    --GMMATmodelFile=${in_gmat} \
#    --varianceRatioFile=${in_var} \
#    --SAIGEOutputFile=${prefix_chr} \
#    --numLinesOutput=2 \
#    --IsOutputAFinCaseCtrl=TRUE \
#    --IsOutputNinCaseCtrl=TRUE \
#    --IsOutputHetHomCountsinCaseCtrl=TRUE \
#    --LOCO=FALSE \
#    ${1:+--condition "$1"}
#  set +x
#}

spa_chromosome_loop() {
   printf 'Markers %s\n' "$1"
   for chr in {7..8}; do
      local prefix_chr="${prefix_iter}_chr${chr}"
      set -x
      Rscript "${step2_SPAtests}"  \
        --vcfFile=${vcf} \
        --vcfFileIndex=${csi} \
        --vcfField="GT" \
        --chrom="chr${chr}" \
        --minMAF=0.0000001 \
        --minMAC=1 \
        --GMMATmodelFile=${in_gmat} \
        --varianceRatioFile=${in_var} \
        --SAIGEOutputFile=${prefix_chr} \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE \
        --IsOutputNinCaseCtrl=TRUE \
        --IsOutputHetHomCountsinCaseCtrl=TRUE \
        --LOCO=FALSE \
        ${1:+--condition "$1"}
      set +x
   done 
}

get_sig_marker() {
  echo $( cat "${prefix_iter}.txt" | \
    awk -v P="${P_cutoff}" '$13 < P' | \
    sort -k 13,13 | \
    cut -d" " -f3 | \
    head -n1)
}




i=0
max=2
marker_list=""
set_up_RSAIGE
while [ $i -lt $max ]
do
    echo "iteration: $i"
    true $(( i++ ))
    
    prefix_iter="${out_prefix}_i${i}"
    spa_chromosome_loop ${marker_list}
    #cat "${prefix_iter}_chr"* > "${prefix_iter}.txt"
    #rm  "${prefix_iter}_chr"*
    marker=$(get_sig_marker)

    echo "echo ${marker}"
    if [ ! -z ${marker} ]; then
        if [ -z ${marker_list} ]; then
            markers=${marker}
        else
            marker_list+=",${marker}"        
        fi        
    else
        echo 'ending loop..'
        ${marker_list} > ${out}
        break    
    fi    
    echo "markers: ${marker_list}"

done



