#!/usr/bin/env bash
#
#$ -N wes_naive_ser
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_naive_ser.log
#$ -e logs/wes_naive_ser.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qe
#$ -t 20
#$ -V

set -o errexit
set -o nounset

module load BCFtools/1.12-GCC-10.3.0

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly in_dir="data/phased/wes/naive"
readonly in_prefix="${in_dir}/ukb_eur_wes_prefilter_200k_naive_phasing_chr${chr}"

files=""
for f in ${in_prefix}.vcf.gz; do 
  make_tabix ${f} "tbi"
  files="${files} ${f}"
done


readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"

readonly out_dir="data/phased/naive"
readonly out_prefix="${out_dir}/ukb_eur_wes_prefilter_200k_naive_phasing_chr${chr}"
readonly out="${out_prefix}.vcf.gz"
readonly trio="${out_prefix}.trio"

mkdir -p ${out_dir}

# count files
n=$( ls -l ${in_prefix}*.vcf.gz | wc -l )

if [ ${n} -gt 1 ] && [ ! -f ${out} ]; then 
  bcftools concat --ligate ${files} -O z -o ${out}
else
  ln -s ${in_prefix}*.vcf.gz ${out}
fi


if [ ! -f ${trio} ]; then
    bcftools +trio-switch-rate ${out} -- -p ${pedigree} > ${trio}
    switch_errors_by_site ${out} ${pedigree}
fi


