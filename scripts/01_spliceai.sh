#!/usr/bin/env bash
#
#$ -N spliceai
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/spliceai.log
#$ -e logs/spliceai.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q short.qc
#$ -t 21
#$ -V

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

#readonly in_dir="data/vep/full"
readonly in_dir="data/unphased/wes_union_calls"
readonly out_dir="data/vep/spliceai"

readonly grch38="/well/lindgren/flassen/software/SpliceAI/spliceai/annotations/grch38.txt"
readonly genome="/well/lindgren/flassen/software/SpliceAI/fasta/hg38.fa.gz"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly tmp_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_rows_chr${chr}"
readonly tmp_file="${tmp_prefix}.vcf.bgz"
readonly out_file="${out_dir}/ukb_wes_200k_full_spliceai_chr${chr}.vcf"

readonly hail_script="scripts/01_spliceai.py"

mkdir -p ${out_dir}

SECONDS=0
if [ ! -f ${tmp_file} ]; then
set_up_hail
set_up_pythonpath_legacy
python3 "${hail_script}" \
     --input_path "${in_file}" \
     --input_type "vcf" \
     --out_prefix "${tmp_prefix}" \
     && print_update "Finished writing rows for chr${chr}" ${SECONDS} \
     || raise_error "Hail writing rows for chr${chr} failed"
fi 

module purge
module load BCFtools/1.12-GCC-10.3.0
make_tabix ${tmp_file} "tbi"

module purge
set_up_tensorflow
spliceai -I ${tmp_file} -O ${out_file} -R ${genome} -A ${grch38}





