#!/usr/bin/env bash
#
#$ -N to_vcf
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/to_vcf.log
#$ -e logs/to_vcf.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc@@short.hga
#$ -t 1-21

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/03_to_vcf.py"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 

readonly in_dir="data/mt/annotated"
readonly input_prefix="${in_dir}/ukb_eur_wes_200k_annot_chr${chr}.mt"
readonly input_type="mt"

readonly out_dir="data/mt/annotated"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_annot_chr${chr}"
readonly out_type="vcf"

SECONDS=0
set_up_hail
set_up_pythonpath_legacy  
set -x
python3 "${hail_script}" \
   --in_file ${input_prefix}\
   --in_type "${input_type}" \
   --out_prefix ${out_prefix} \
   --out_type "${out_type}" \
   && print_update "Finished annotating MatrixTables chr${chr}" ${SECONDS} \
   || raise_error "Annotating MatrixTables for chr${chr} failed"
set +x 


if [ ${out_type} == "vcf" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix}.vcf.bgz" "csi"
fi






