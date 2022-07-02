#!/usr/bin/env bash
#
#$ -N mac_by_phenotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/mac_by_phenotypes.log
#$ -e logs/mac_by_phenotypes.errors.log
#$ -P lindgren.prjc
#$ -q short.qc
#$ -pe shmem 10
#$ -t 20-22
#$ -V


set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly rscript="scripts/conditional/rare/02_mac_by_phenotypes.R"

readonly chr="${SGE_TASK_ID}"
readonly pheno_dir="data/phenotypes"
readonly in_dir="data/conditional/rare/combined"
readonly out_dir="data/conditional/rare/combined"

readonly pheno_path="${pheno_dir}/curated_covar_phenotypes_cts_200k.tsv"
readonly input_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_misseense_AC_cts"

readonly phenotypes=$(cat "${pheno_dir}/filtered_phenotypes_cts_manual.tsv" | tr "\n" "," | sed 's/\(.*\),/\1 /' )

set_up_rpy
Rscript "${rscript}" \
   --in_vcf ${input_path} \
   --out_prefix ${out_prefix} \
   --subset_phenotypes ${phenotypes} \
   --phenotypes ${pheno_path}



