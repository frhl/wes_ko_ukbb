#!/usr/bin/env bash
#
#$ -N index
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/index.log
#$ -e logs/index.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qe
#$ -t 1-22

source utils/vcf_utils.sh

module load BCFtools/1.12-GCC-10.3.0

# paths (ptvs_damaging_missense/)
readonly chr=${SGE_TASK_ID}
readonly in_dir="derived/knockouts/all/ptvs_damaging_missense/"
readonly in_file="${in_dir}/ukb_wes_200k_maf002_miss005_ptv_chr${chr}_ko.vcf.bgz"

# index
make_tabix "${in_file}" "csi"





