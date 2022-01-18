#!/usr/bin/env bash
#
#$ -N single_predictor
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/single_predictor.log
#$ -e logs/single_predictor.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qa@@short.hga
#$ -t 21
#$ -V


readonly bed_dir="data/prs/imp/"
readonly out_dir="data/prs/single_predictor"
readonly grm_dir="data"
readonly pheno_dir="data/phenotypes"
readonly covar_dir="data/phenotypes"

readonly ldak_dir="/well/lindgren/flassen/software/ldak"
readonly ldak="${ldak_dir}/ldak5.2.linux"

readonly chr=$( get_chr ${SGE_TASK_ID} )
readonly bed="${out_dir}/"
readonly grm="${grm_dir}/"
readonly phenos="${pheno_dir}/curated_phenotypes.tsv"
readonly covars="${covar_dir}/curated_phenotypes.tsv"
readonly index="1"

readonly out_prefix="${out_dir}/ukb_imp500k_single_predictor${chr}"

readonly sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

mkdir -p ${out_dir}

./ "${ldak}" \
  --bfile ${bed} \
  --grm ${grm} \
  --ch4 ${chr} \
  --pheno ${phenos} \
  --mpheno ${index} \
  --linear \

# what is the outfile argument?





