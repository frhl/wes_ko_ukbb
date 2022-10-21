#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=annotate
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/annotate.log
#SBATCH --error=logs/annotate.errors.log
#SBATCH --partition=long
#SBATCH --cpus-per-task 5
#SBATCH --array=21
#SBATCH --requeue

# Note: long.qc@@long.hga with 4 slots required to run full pipeline
# -q long.qc@@long.hga

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/03_annotate.py"

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} ) 
readonly in_dir="data/mt/union"
readonly input_prefix="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.mt"

readonly out_dir="data/mt/annotated/test"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly out_type="vcf"
readonly out="${out_prefix}.vcf.bgz"
readonly out_trio="${out_prefix}.trio"

readonly annotation_table="data/vep/hail/ukb_wes_200k_chr${chr}_vep.ht"

readonly pedigree_dir="/well/lindgren/UKBIOBANK/nbaya/resources"
readonly pedigree="${pedigree_dir}/ukb11867_pedigree.fam"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly final_variant_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

mkdir -p ${out_dir}



if [ ! -f ${out} ]; then
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
     --out_type "${out_type}" \
     --dbsnp_path "155" \
     --annotate_rsid \
     --annotate_snp_id \
     && print_update "Finished annotating MatrixTables chr${chr}" ${SECONDS} \
     || raise_error "Annotating MatrixTables for chr${chr} failed"
  set +x 
fi





if [ -f "${out}" ]; then
  module purge
  module load BCFtools/1.12-GCC-10.3.0
  make_tabix "${out_prefix}.vcf.bgz" "csi"
fi

make_tabix "${out_prefix}.vcf.bgz" "tbi"



