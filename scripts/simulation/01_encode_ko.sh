#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=encode_ko
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/encode_ko.log
#SBATCH --error=logs/encode_ko.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 3
#SBATCH --array=22
#
#$ -N encode_ko
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/encode_ko.log
#$ -e logs/encode_ko.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qa
#$ -t 22
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/simulation/01_encode_ko.py"
readonly spark_dir="data/tmp/spark_dir"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly n_samples="100k"
readonly in_dir="data/simulation/mt"
readonly in_file="${in_dir}/ukb_eur_${n_samples}_chr${chr}.mt"
readonly in_type="mt"

readonly out_dir="data/simulation/knockouts"
readonly out_prefix="${out_dir}/ukb_eur_${n_samples}_encoded_norm_knockouts_chr${chr}"
readonly out_type="vcf"

# Allele frequency thresholds to filter on
readonly maf_min="0"
readonly maf_max="5e-2"

# Discard probabilistic encodings? I.e. unphased singletons
#readonly discard_prob_dosages="Y"

# What group of SNPs should be kept
readonly in_category="pLoF,damaging_missense"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

# Sample WES
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
      --chrom ${chr} \
      --in_prefix ${in_file} \
      --in_type ${in_type} \
      --csqs_category ${in_category} \
      --use_loftee \
      --maf_max ${maf_max} \
      --maf_min ${maf_min} \
      --out_prefix ${out_prefix} \
      --out_type ${out_type} 




