#!/usr/bin/env bash
#
#$ -N prefilter_phenotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter_phenotypes.log
#$ -e logs/prefilter_phenotypes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q short.qc
#$ -t 21
#$ -tc 10
#$ -V

set -o errexit
set -o nounset

module purge
source utils/bash_utils.sh

readonly rscript="scripts/conditional/rare/02_prefilter_phenotypes.R"

readonly chr="${SGE_TASK_ID}"
readonly pheno_dir="data/phenotypes"
readonly in_dir="data/conditional/rare/combined"
readonly out_dir="data/conditional/rare/combined"

readonly pheno_cts_path="${pheno_dir}/curated_covar_phenotypes_cts_200k.tsv"

readonly in_vcf="${in_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly tmp_vcf="${in_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.txt"
readonly tmp_vcf_gz="${tmp_vcf}.gz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense_test"
readonly out_mrg="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense_test"
readonly covar_path="${pheno_dir}/covars1.csv"

readonly phenotypes_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
readonly phenotypes_bin="${pheno_dir}/filtered_phenotypes_binary_header.tsv"

mkdir -p ${out_dir}

# need to generate temporary uncompressed VCF
# otherwise, fread will have to read chunks from 
# zcat which is ineffecient (?)
if [ ! -f "${tmp_vcf_gz}" ]; then
  >&2 echo "checking $in_vcf"
  zcat "${in_vcf}" | grep -v "##" > ${tmp_vcf}
  echo -e "\n" >> ${tmp_vcf}
  gzip ${tmp_vcf}
  rm -f ${tmp_vcf}
fi

submit_binary(){
  local trait="binary"
  local pheno_list="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  local pheno_file="${pheno_dir}/curated_covar_phenotypes_binary_200k.tsv" 
  submit_qc_job ${pheno_list} ${pheno_file} ${trait}
}

submit_qc_job() {
  local pheno_list=${1}
  local pheno_file=${2}
  local trait=${3}
  local qsub_main="_pref_c${chr}_${trait}"
  local qsub_mrg="_mrg_c${chr}"
  local pheno_list_csv=$(cat ${pheno_list} | tr "\n" ",")

  set_up_rpy
  Rscript "${rscript}" \
    --phenotypes ${pheno_list_csv} \
    --pheno_file ${pheno_file} \
    --in_vcf ${tmp_vcf_gz} \
    --covariates ${covar_path} \
    --out_prefix ${out_prefix}

}

submit_binary


