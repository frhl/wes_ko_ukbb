#!/usr/bin/env bash
#
#$ -N annotate_gnomad_vep_qc
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/annotate_gnomad_vep_qc.log
#$ -e logs/annotate_gnomad_vep_qc.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc@@short.hge
#$ -t 1-24

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories

readonly spark_dir="data/tmp/spark"
readonly out_dir="/well/lindgren/flassen/ressources/gnomad/gnomad_v2_liftover/hail_vep"
readonly in_dir="/well/lindgren/flassen/ressources/gnomad/gnomad_v2_liftover/exomes"

# hail script
readonly hail_script="scripts/QC/07_annotate_gnomad_vep.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/gnomad.exomes.r2.1.1.sites.${chr}.liftover_grch38.vcf.bgz"

# output path
readonly out_prefix="${out_dir}/gnomad.exomes.r2.1.1.sites.${chr}.liftover_grch38"

mkdir -p ${out_dir}
set_up_hail
set_up_vep
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "vcf" \
     --out_prefix ${out_prefix}

