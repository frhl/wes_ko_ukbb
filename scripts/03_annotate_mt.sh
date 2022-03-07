#!/usr/bin/env bash
#
#$ -N annotate_mt
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/annotate_mt.log
#$ -e logs/annotate_mt.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q long.qc@@long.hga
#$ -t 14-15

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/03_annotate_mt.py"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_dir="data/mt/union"
readonly input_prefix="${in_dir}/ukb_eur_wes_200k_union_chr${chr}.mt"

readonly out_dir="data/mt/annotated"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_annot_chr${chr}"
readonly out="${out_prefix}.mt"

readonly annotation_table="data/vep/hail/ukb_wes_200k_chr${chr}_vep.ht"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

mkdir -p ${out_dir}

if [ ! -f "${out}/_SUCCESS" ]; then
  SECONDS=0
  set_up_hail
  set_up_pythonpath_legacy  
  set -x
  python3 "${hail_script}" \
     --in_file ${input_prefix}\
     --in_type "mt" \
     --input_annotation_path ${annotation_table}\
     --final_sample_list ${final_sample_list} \
     --final_variant_list ${final_variant_list}\
     --out_prefix ${out_prefix} \
     --out_type "mt" \
     --dbsnp_path "155" \
     --annotate_rsid \
     --annotate_snp_id \
     && print_update "Finished annotating MatrixTables chr${chr}" ${SECONDS} \
     || raise_error "Annotating MatrixTables for chr${chr} failed"
  set +x 
else
  print_update "file ${out} already exists. Skipping!"
fi






