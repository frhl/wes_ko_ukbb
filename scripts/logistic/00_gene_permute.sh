#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=gene_permute
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/gene_permute.log
#SBATCH --error=logs/gene_permute.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22


set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly rscript="scripts/logistic/00_gene_permute.R"

readonly chr=$( get_array_task_id )
readonly seed=3

readonly input_dir="data/knockouts/alt/pp90/combined"
readonly input_path="${input_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense.vcf.bgz"

readonly out_dir="data/knockouts/alt/pp90/combined/shuffled/seed${seed}"
readonly out_prefix="${out_dir}/shuffled_ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense"

mkdir -p ${out_dir}

# create shuffled VCF
set_up_rpy
Rscript ${rscript} \
  --chrom "chr${chr}" \
  --input_path ${input_path} \
  --out_prefix ${out_prefix} \
  --seed ${seed} \

# zip file and tabix  
module purge
module load BCFtools/1.12-GCC-10.3.0
bgzip "${out_prefix}.vcf"
rm -f "${out_prefix}.vcf"
make_tabix "${out_prefix}.vcf.gz" "csi"




