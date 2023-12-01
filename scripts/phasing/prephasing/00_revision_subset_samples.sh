#!/usr/bin/env bash
#
# @description create matrix table of random samples to be used for read-backed phasing
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=revision_subset_samples
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/revision_subset_samples.log
#SBATCH --error=logs/revision_subset_samples.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/vcf_utils.sh

readonly rscript="scripts/phasing/prephasing/00_revision_subset_samples.R"

readonly input_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/readbacked/data/phased/subset/50k"
readonly input_path="${input_dir}/UKB.wes.200k.chr21.50k.b1of4.vcf.gz"

readonly qced_samples_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/samples"
readonly qced_samples="${qced_samples_dir}/UKB.chr21.samples.eur.txt"

readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"

readonly n_samples=10000
readonly seed=1462

readonly out_dir="data/prephased/wes_union_calls/revision/10k"
readonly out_prefix="${out_dir}/wes.revision.10k.${seed}"

mkdir -p ${out_dir}

# get overlapping samples
module load BCFtools
readonly vcf_samples="${out_prefix}.samples_in_vcf"

bcftools query -l ${input_path} > ${vcf_samples}

echo "wc -l vcf_samples"
wc -l ${vcf_samples}

# load R and get subset
module purge && module load R
Rscript ${rscript} \
  --sample_file_qc ${qced_samples} \
  --sample_file_vcf ${vcf_samples} \
  --pedigree_file ${pedigree} \
  --n_samples ${n_samples} \
  --out_prefix ${out_prefix} \
  --seed ${seed}

# clean up
rm ${vcf_samples}

