#!/usr/bin/env bash
#
#$ -N prefilter
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prefilter.log
#$ -e logs/prefilter.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc@@short.hga
#$ -t 21


module load BCFtools/1.12-GCC-10.3.0

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly in_dir="data/unphased/post-qc/"
readonly out_dir="data/phased/test-phasing"

readonly chr="${SGE_TASK_ID}"
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
readonly out_prefix="${out_dir}/ukb_wes_200k_filtered_maf001_chr${chr}"
readonly hail_script="scripts/debug/10_prefilter.py"

set_up_hail
set_up_pythonpath_legacy
mkdir -p ${out_dir}

python3 "${hail_script}" \
   --input_path ${in_file} \
   --input_type "vcf" \
   --min_maf 0.01 \
   --out_prefix ${out_prefix} 
print_update "Hail finished writing."
make_tabix "${out_prefix}.vcf.bgz" "tbi"





