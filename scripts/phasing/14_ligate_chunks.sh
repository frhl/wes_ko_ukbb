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
#SBATCH --array=21

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly rscript="scripts/phasing/_sort_chunks.R"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} )
readonly in_dir="data/phased/wes_scaffold_calls/200k_from_500k/trimmed"
readonly in_prefix="${in_dir}/ukb_wes_scaffold_calls_200k_from_500k_chr${chr}_trim.1of1"
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

readonly out_dir="data/phased/wes_scaffold_calls/200k_from_500k/ligated"
readonly out_prefix="${out_dir}/ukb_wes_scaffold_calls_200k_fromm_500k_chr${chr}"
readonly out="${out_prefix}.vcf.bgz"

mkdir -p ${out_dir}


if [ ! -f "${out}" ]; then
  if [ ${n} -gt 1 ]; then 
    echo "Ligating ${n} files for chromosome ${chr}.."
    bcftools concat --ligate ${files} -O z -o ${out}
  else
    ln -s "${PWD}/${in_prefix}"*.vcf.bgz ${out}
  fi
fi

if [ ! -f "${out}.tbi" ]; then
  make_tabix ${out}
fi




