#!/usr/bin/env bash
#
#$ -N wes_union_call_ser
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_union_call_ser.log
#$ -e logs/wes_union_call_ser.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 20-22
#$ -V


readonly in_dir="data/phased/wes/union_call"
readonly out_dir="data/switch"
readonly fam_dir="/well/lindgren/UKBIOBANK/nbaya/resources"

readonly chr="${SGE_TASK_ID}"
readonly in_phased="${in_dir}/ukb_wes_union_calls_200k_chr${chr}.vcf.gz"
readonly fam_file="${fam_dir}/ukb11867_pedigree.fam"

readonly var_file="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.var"
readonly ser_file="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.ser"
readonly out_file="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.txt"

mkdir -p ${out_dir}

switch_errors_by_site() {

  module load R
  module load BCFtools/1.12-GCC-10.3.0
  export BCFTOOLS_PLUGINS="/well/lindgren/flassen/software/bcf/bcftools-1.12/plugins"

  local vcf=${1}
  local pedigree=${2}
  local rscript="utils/aggr_ser.R" 
  local ser="${vcf%.*}.ser"
  local var="${vcf%.*}.var"
  local out="${vcf%.*}.txt"
  set -x
  bcftools +trio-switch-rate "${vcf}" -- -p "${pedigree}" -c 0 > "${ser_file}"
  bcftools +fill-tags "${vcf}" | \
  bcftools query -f '%ID %CHROM %POS %REF %ALT %MAF %AF %AC %AN %HWE\n'  -o "${var_file}"
  Rscript ${rscript} \
    --switch-error-file ${ser_file} \
    --variant-file ${var_file} \
    --outfile ${out_file}
  set +x
  rm ${ser} ${var}

}
