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
#SBATCH --array=2-22
#
#$ -N ligate_chunks
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/ligate_chunks.log
#$ -e logs/ligate_chunks.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qc
#$ -t 20-22
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/vcf_utils.sh

readonly rscript="scripts/phasing/phasing/_sort_chunks.R"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly PWD=$(pwd)

# eagle2
#readonly in_dir="data/phased/wes_union_calls/200k/eagle2/trimmed"
#readonly in_prefix="${in_dir}/ukb_wes_union_calls_eagle2_200k_chr${chr}"
#readonly in_trim="${in_prefix}_trim_trim"
#readonly out_dir="data/phased/wes_union_calls/200k/eagle2/ligated"
#readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
#readonly out="${out_prefix}.vcf.bgz"

# SHAPEIT4
readonly in_dir="data/phased/wes_union_calls/200k/shapeit4/trimmed"
readonly in_prefix="${in_dir}/ukb_wes_union_calls_shapeit4_200k_chr${chr}"
readonly in_trim="${in_prefix}_trim_trim"
readonly out_dir="data/phased/wes_union_calls/200k/shapeit4/ligated"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly out="${out_prefix}.vcf.bgz"

# SHAPEIT5
#readonly in_dir="data/phased/wes_union_calls/200k/shapeit5/trimmed"
#readonly in_prefix="${in_dir}/ukb_wes_union_calls_shapeit5_200k_chr${chr}"
#readonly in_trim="${in_prefix}_trim_trim"
#readonly out_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
#readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"
#readonly out="${out_prefix}.vcf.bgz"


set_up_rpy
module load BCFtools/1.12-GCC-10.3.0
for f in ${in_prefix}*.vcf.bgz; do
  make_tabix ${f} "tbi"
done

readonly files=$( Rscript ${rscript} --in_prefix ${in_trim})
readonly n=$(echo ${files} | tr " " "\n" | wc -l)

echo $files
echo "\nNote: Chunks found ${n}"

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




