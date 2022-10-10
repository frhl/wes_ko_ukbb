#!/usr/bin/env bash
#
# @description combine chunks in a phase-aware manner into full chromosomes.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=ligate_chunks
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/ligate_chunks.log
#SBATCH --error=logs/ligate_chunks.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=22

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh

set -x

readonly rscript="scripts/phasing/_sort_chunks.R"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly in_dir="data/phased/wes_union_calls/trimmed"
readonly in_prefix="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly in_trim="${in_prefix}_trim"

set_up_rpy
module load BCFtools/1.12-GCC-10.3.0
for f in ${in_prefix}*.vcf.bgz; do 
  make_tabix ${f} "tbi"
done

readonly files=$( Rscript ${rscript} --in_prefix ${in_trim})
readonly n=$(echo ${files} | tr " " "\n" | wc -l)

echo $files
echo "\nNote: Chunks found ${n}"

readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"

readonly out_dir="data/phased/wes_union_calls/ligated"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly out="${out_prefix}.vcf.bgz"
readonly trio="${out_prefix}.trio"

mkdir -p ${out_dir}


if [ ! -f ${out} ]; then
  if [ ${n} -gt 1 ]; then 
    bcftools concat --ligate ${files} -O z -o ${out}
  else
    ln -s "${PWD}/${in_prefix}"*.vcf.bgz ${out}
  fi
fi

if [ ! -f ${trio} ]; then
    bcftools +trio-switch-rate ${out} -- -p ${pedigree} > ${trio}
    switch_errors_by_site ${out} ${pedigree}
fi


