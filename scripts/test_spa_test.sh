#!/usr/bin/env bash
#
#
#$ -N test_spa_test
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/test_spa_test.log
#$ -e logs/test_spa_test.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly phenotype="Alanine_aminotransferase_residual"
#in_vcf="data/knockouts/alt/ukb_eur_wes_200k_chr1_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
#in_csi="data/knockouts/alt/ukb_eur_wes_200k_chr1_maf0to5e-2_pLoF_damaging_missense.vcf.bgz.csi"

readonly in_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/permute/permutations/chrs/chr1/ENSG00000162654"
#readonly in_vcf="${in_dir}/ukb_eur_wes_200k_chr1_maf0to5e-2_pLoF_damaging_missense_two_markers.vcf.gz"
readonly in_vcf="${in_dir}/ukb_eur_wes_200k_chr1_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly in_csi="${in_vcf}.csi"
#readonly in_vcf="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chr1_ENSG00000162654_1.vcf.gz"
#readonly in_csi="${in_vcf}.csi"
readonly in_gmat="data/saige/output/cts/step1/ukb_wes_200k_Alanine_aminotransferase_residual.rda"
readonly in_var="data/saige/output/cts/step1/ukb_wes_200k_Alanine_aminotransferase_residual.varianceRatio.txt"
readonly min_mac=4
readonly out_prefix="${in_dir}/test_new_saige"
readonly chr=1
readonly in_markers=""


# chromosome specified at in_gmat and in_var
# when off chromosome PRS will beused
readonly gmat=$(echo ${in_gmat} | sed -e "s/CHR/${chr}/g")
readonly var=$(echo ${in_var} | sed -e "s/CHR/${chr}/g")

>&2 echo "${gmat} and ${var} with ${phenotype}"

readonly var_bytes=$( file_size ${var} )
readonly gmat_bytes=$( file_size ${gmat} )

readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly markers=$(cat $(echo ${in_markers} | sed -e "s/CHR/${chr}/g"))
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")

#readonly threads=$(( ${NSLOTS}-1 ))
readonly threads=1
readonly step2_SPAtests="utils/saige/step2_SPAtests.R"

spa_test() {
  echo "var_bytes=${var_bytes} at ${var}"
  echo "gmat_bytes=${gmat_bytes} at ${gmat}"
  if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then 
  SECONDS=0
  Rscript "${step2_SPAtests}" \
     --vcfFile=${vcf} \
     --vcfFileIndex=${csi} \
     --vcfField="DS" \
     --chrom="chr${chr}" \
     --minMAF=0.0000001 \
     --minMAC=${min_mac} \
     --GMMATmodelFile=${gmat} \
     --varianceRatioFile=${var} \
     --SAIGEOutputFile=${out} \
     --LOCO=FALSE\
     && print_update "Finished saddle-point approximation. Writing to ${out}" ${SECONDS} \
     || raise_error "Saddle-point approximation for chr${chr} failed"
  else
    raise_error "${var} or ${gmat} does not contain any bytes!"
  fi
}

if [ ! -f ${out} ]; then
   set_up_RSAIGE
   spa_test
else
  >&2 echo "${out} already exists. Skipping.."
fi 


