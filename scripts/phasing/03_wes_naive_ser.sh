#!/usr/bin/env bash
#
#$ -N wes_naive_ser
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_naive_ser.log
#$ -e logs/wes_naive_ser.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 20
#$ -V

module load BCFtools/1.12-GCC-10.3.0

readonly in_dir="data/phased/wes/naive"
readonly out_dir="data/switch"
readonly fam_dir="/well/lindgren/UKBIOBANK/nbaya/resources"

readonly chr="${SGE_TASK_ID}"
readonly in_phased="${in_dir}/ukb_wes_200k_naive_phasing_chr${chr}.vcf.gz"
readonly fam_file="${fam_dir}/ukb11867_pedigree.fam"
readonly rscript="scripts/evaluation/aggregate_switch_errors.R"

readonly prefix="${out_dir}/ukb_wes_200k_naive_phasing_chr${chr}"
readonly var_file="${prefix}.var"
readonly ser_file="${prefix}.ser"
readonly out_file="${prefix}.txt"

export BCFTOOLS_PLUGINS="/well/lindgren/flassen/software/bcf/bcftools-1.12/plugins"

mkdir -p ${out_dir}

if [ ! -f "${ser_file}" ]; then
  set -x
  bcftools +trio-switch-rate "${in_phased}" -- -p "${fam_file}" -c 0 > "${ser_file}"
  set +x
fi 

if [ ! -f "${var_file}" ]; then
  set -x
  bcftools +fill-tags "${in_phased}" | \
  bcftools query -f '%ID %CHROM %POS %REF %ALT %MAF %AF %AC %AN %HWE\n'  -o "${var_file}"
  set +x
fi

module purge
module load R
Rscript ${rscript} \
  --switch-error-file ${ser_file} \
  --variant-file ${var_file} \
  --outfile ${out_file}



