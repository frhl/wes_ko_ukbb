#!/usr/bin/env bash
#
#$ -N mac_by_phenotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/mac_by_phenotypes.log
#$ -e logs/mac_by_phenotypes.errors.log
#$ -P lindgren.prjc
#$ -q short.qc
#$ -pe shmem 4
#$ -t 21
#$ -V


set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/rare/02_mac_by_phenotypes.py"

readonly chr="${SGE_TASK_ID}"
readonly pheno_dir="data/phenotypes"
readonly in_dir="data/conditional/rare/combined/test"
readonly out_dir="data/conditional/rare/combined/test"

readonly pheno_path="${pheno_dir}/curated_covar_phenotypes_cts_200k.tsv"
readonly input_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_misseense_by_cts_pheno"
readonly input_type="vcf"
readonly out_type="vcf"

readonly phenotypes=$(cat "${pheno_dir}/filtered_phenotypes_cts_manual.tsv" | tr "\n" "," | sed 's/\(.*\),/\1 /' )

set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --input_path ${input_path} \
   --input_type ${input_type} \
   --out_prefix ${out_prefix} \
   --out_type ${out_type} \
   --pheno_path ${pheno_path} \
   --phenotypes ${phenotypes}



