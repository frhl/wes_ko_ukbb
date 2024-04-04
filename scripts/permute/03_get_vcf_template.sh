#!/usr/bin/env bash
#
# @description for each gene across all phenotypes figure out how many permutations are needed
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_vcf_template
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/get_vcf_template.log
#SBATCH --error=logs/get_vcf_template.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

source utils/vcf_utils.sh
module load BCFtools/1.12-GCC-10.3.0

# need any VCF with knockouts to get order of samples
readonly in_dir="data/knockouts/alt/pp90/combined"
readonly in_vcf="${in_dir}/ukb_eur_wes_200k_chr21_pLoF.vcf.bgz"

readonly out_dir="data/permute/overview"
readonly out_file="${out_dir}/sample_order.txt"

# get sample order that is used to create
mkdir -p ${out_dir}


# NOTE: there is an error in SLURM/BCFTOOLS causing this
# to fail when submitted as a job. You might need to do it manuall (Mar-2023)
# update: i think this issue stems from using 'epyc' nodes
echo "$( bcftools query -l ${in_vcf} )" > ${out_file}



