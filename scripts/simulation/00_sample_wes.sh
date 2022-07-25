#!/usr/bin/env bash
#
#$ -N sample_wes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/sample_wes.log
#$ -e logs/sample_wes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 3
#$ -q short.qa
#$ -t 21
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly hail_script="scripts/simulation/00_sample_wes.py"
readonly spark_dir="data/tmp/spark_dir"

readonly chr=$( get_chr ${SGE_TASK_ID} )

readonly in_dir="data/mt/annotated"
readonly in_file="${in_dir}/ukb_eur_wes_200k_annot_chr${chr}.mt"
readonly in_type="mt"

readonly ko_dir="data/knockouts/alt"
readonly ko_file="${ko_dir}/ukb_eur_wes_200k_chr${chr}_maf0to5e-2_pLoF_damaging_missense.vcf.bgz"
readonly ko_type="vcf"

readonly out_dir="data/simulation/mt"
readonly out_prefix="${out_dir}/ukb_eur_100k_chr${chr}"
readonly out_type="vcf"

readonly seed="1995"

mkdir -p ${spark_dir}
mkdir -p ${out_dir}

# Sample WES
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
   --in_prefix "${in_file}"\
   --in_type "mt" \
   --random_samples 100000 \
   --random_seed 1995 \
   --filter_to_unrelated_using_kinship_coef \
   --out_prefix "${out_prefix}" \
   --out_type "${out_type}" 

# Sample corresponding KO file




