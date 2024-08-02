#!/usr/bin/env bash
#
# @description generate files of genotyped calls
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=split_calls
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/split_calls.log
#SBATCH --error=logs/split_calls.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 3
#SBATCH --array=20
#
#$ -N split_calls
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/split_calls.log
#$ -e logs/split_calls.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qe
#$ -t 1
#$ -V


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/phasing/prefilter/split_parents.py"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly tranche="200k"

readonly in_dir="data/unphased/calls/prefilter/${tranche}"
readonly in_file="${in_dir}/ukb_prefilter_calls_${tranche}_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/unphased/calls/prefilter/${tranche}"
readonly out_prefix="${out_dir}/ukb_split_calls_${tranche}_chr${chr}"
readonly out_prefix_parents="${out_prefix}_parents"
readonly out_prefix_no_parents="${out_prefix}_no_parents"
readonly out_type="vcf"

# samples overlapping exomes and genotypes
readonly parents_dir="data/unphased/overlap"
readonly parents_path="${parents_dir}/ukb_calls_wes_samples_parents.txt"

# fasta reference
readonly fasta="/well/lindgren/flassen/ressources/genome_reference/broad/Homo_sapiens_assembly38.fasta"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

if [ ! -f "${out_prefix_no_parents}.vcf.bgz" ]; then
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
     --input_path ${in_file} \
     --input_type ${in_type} \
     --parents_path ${parents_path} \
     --out_prefix "${out_prefix}" \
     --out_type "${out_type}" 
fi

if [ ! -f "${out_prefix_no_parents}.vcf.bgz.tbi" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix_no_parents}.vcf.bgz" "tbi"
fi

if [ -f "${out_prefix_parents}.vcf.bgz" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix_parents}.vcf.bgz" "tbi"
  bcftools +fixref "${out_prefix_no_parents}.vcf.bgz" -- -f ${fasta} > "${out_prefix_no_parents}.fasta.txt"
fi



