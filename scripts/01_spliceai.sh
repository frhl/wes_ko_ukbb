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
source utils/bash_utils.sh

readonly in_dir="data/unphased/wes_union_calls"
readonly out_dir="data/vep/spliceai"

readonly genome="/well/indgren/flassen/software/SpliceAI/spliceai/annotations/grch38.txt"
readonly grch38="/well/lindgren/flassen/software/SpliceAI/fasta/hg38.fa.gz"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly infile="${in_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.bgz"
readonly outfile="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}.vcf.gz"

set_up_spliceai
spliceai -I ${infile} -O ${outfile} -R ${genome} -A ${grch38}




