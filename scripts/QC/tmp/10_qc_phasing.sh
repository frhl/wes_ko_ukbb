#!/usr/bin/env bash
#
#$ -N unphased_hardcalls
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/create_hardcalls_unphased.log
#$ -e logs/create_hardcalls_unphased.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 8
#$ -q short.qc@@short.hge
#$ -t 1-24

source utils/qsub_utils.sh
source utils/hail_utils.sh

# directories
readonly in_dir="data/mt"
readonly ref_dir="/well/1000G/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

# hail script
readonly hail_script="scripts/QC/10_qc_phasing.py"

# input path
readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly in_file="${in_dir}/ukb_wes_200k_filtered_chr${chr}.mt"
readonly in_ref="${ref_dir}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"
readonly in_phased="${in_dir}/ukb_wes_200k_annotated_chr${chr}.mt"

# output path
readonly out_prefix="${out_dir}/ukb_wes_200k_chr${chr}"

mkdir -p ${out_dir}
set_up_hail
set_up_pythonpath_legacy  
python3 "${hail_script}" \
     --input_path ${in_file}\
     --input_type "mt" \
     --input_reference ${in_ref} \
     --out_prefix ${out_prefix}


