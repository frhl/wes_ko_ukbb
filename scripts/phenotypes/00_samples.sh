#!/usr/bin/env bash

#SBATCH -A lindgren.prj
#SBATCH -J samples
#SBATCH -D /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/samples.log
#SBATCH --error=logs/samples.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#
#
#$ -N samples
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/samples.log
#$ -e logs/samples.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly rscript="scripts/00_samples.R"
readonly hail_script="scripts/00_samples.py"
readonly spark_dir="data/tmp/spark"

readonly vcf_sample_dir="data/phased/wes_union_calls/200k/shapeit5/ligated"
readonly vcf_sample="${vcf_sample_dir}/ukb_wes_union_calls_200k_chr21.vcf.bgz"
readonly qc_palmer="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list"

readonly out_dir="data/phenotypes/samples"
readonly sample_list_qc="${out_dir}/sample_list_final.qc.txt"
readonly sample_list_imputed="${out_dir}/sample_list_imputed.txt"
readonly sample_list_phased="${out_dir}/sample_list_phased.txt"
readonly sample_list_calls="${out_dir}/sample_list_calls.txt"
readonly sample_list_final="${out_dir}/ukb_wes_ko.imputed.qc.samples"
readonly sample_list_final_unrel="${out_dir}/ukb_wes_ko.imputed.qc.unrelated.samples"
readonly sample_list_phased_nfe="${out_dir}/ukb_wes_ko.qc.nfe.samples"

mkdir -p ${out_dir}
mkdir -p ${spark_dir}

# link final list
ln -s ${qc_palmer} ${sample_list_qc}

# get samples that have been phased (without parents)
module load BCFtools/1.12-GCC-10.3.0
bcftools query -l "${vcf_sample}" > "${sample_list_phased}"
module purge

# export genotype array samples
module purge
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --output_path "${sample_list_calls}" \
   --export_genotyped_samples

# export imputed samples
python3 "${hail_script}" \
   --output_path "${sample_list_imputed}" \
   --remove_withdrawn \
   --export_imputed_samples

# export final list of samples (all)
python3 "${hail_script}" \
   --output_path "${sample_list_final}" \
   --extract1 ${sample_list_qc} \
   --extract2 ${sample_list_phased} \
   --remove_withdrawn \
   --export_imputed_samples

# export final list of samples (unrelated)
python3 "${hail_script}" \
   --output_path "${sample_list_final_unrel}" \
   --extract1 ${sample_list_qc} \
   --extract2 ${sample_list_phased} \
   --get_unrelated \
   --remove_withdrawn \
   --export_imputed_samples

# export samples that are NFE and phased (i.e. those
# we end up using for KO analysis)
module purge
set_up_rpy
Rscript ${rscript} \
  --input_path1 ${sample_list_phased} \
  --input_path2 ${sample_list_qc} \
  --outfile ${sample_list_phased_nfe}




