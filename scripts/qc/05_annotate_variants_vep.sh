#!/usr/bin/env bash
#
#$ -N annotate_variants_vep_qc
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/annotate_variants_vep_qc.log
#$ -e logs/annotate_variants_vep_qc.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc@@short.hge
#$ -t 1-24

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories

readonly in_dir="data/unphased/post-qc"
readonly spark_dir="data/tmp/spark"
readonly out_dir="data/variants"
readonly gnomad_dir="/well/lindgren/flassen/ressources/gnomad/gnomad_v2_liftover/exomes"
readonly vep_dir="data/vep/full"

# hail script
readonly hail_script="scripts/QC/05_annotate_variants_vep.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly gnomad="${gnomad_dir}/gnomad.exomes.r2.1.1.sites.${chr}.liftover_grch38.vcf.bgz"
readonly vep="${vep_dir}/ukb_wes_200k_full_vep_chr${chr}.vcf"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_filtered_chr${chr}"

mkdir -p ${out_dir}
set_up_hail
set_up_vep
set_up_pythonpath  
python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "mt" \
     --out_prefix ${out_prefix}\
     --gnomad_path ${gnomad}\
     --vep_path ${vep}


