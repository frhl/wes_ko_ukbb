#!/usr/bin/env bash
#
#$ -N ligate_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ligate_chunks.log
#$ -e logs/ligate_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 5
#$ -q short.qe
#$ -t 1-3
#$ -V

set -o errexit
set -o nounset

module load BCFtools/1.12-GCC-10.3.0

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly main_dir="data/phased/wes_union_calls/chunks/final"
readonly in_dir="${main_dir}/ukb_eur_wes_union_calls_200k_chr${chr}-16xshort.qe"
readonly in_prefix="shapeit4_prs100000_pro25000_mprs100000"

files=""
for f in ${in_dir}/*.vcf.gz; do 
  make_tabix ${f} "tbi"
  files="${files} ${f}"
done

readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"

readonly out_dir="data/phased/wes_union_calls/ligated"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly out="${out_prefix}.vcf.bgz"
readonly trio="${out_prefix}.trio"

mkdir -p ${out_dir}

echo "chr${chr}: ${files}"

if [ ! -f ${out} ]; 
  set -x
  bcftools concat --ligate ${files} -O z -o ${out}
  set +x
fi

if [ ! -f ${trio} ]; then
  bcftools +trio-switch-rate ${out} -- -p ${pedigree} > ${trio}
  switch_errors_by_site ${out} ${pedigree}
fi





