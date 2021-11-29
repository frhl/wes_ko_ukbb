#!/usr/bin/env bash
#
# extract genotypes in regions near genes that are significant in primary analysis
#
#$ -N filter_genotypes
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/filter_genotypes.log
#$ -e logs/filter_genotypes.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc
#$ -t 3

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly pheno_dir="data/phenotypes"
readonly gene_dir="derived/tables/gene_intervals"
readonly out_dir="data/conditional/common"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly gene_table="${gene_dir}/211111_wes200k_saige_wes_saige_sig_genes_intervals.txt"
readonly pheno_list="${pheno_dir}/UKBB_WES200k_binary_phenotypes_header.txt"
readonly hail_script="scripts/conditional/02_filter_genotypes.py"

readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

readonly out_prefix="${out_dir}/211111_ukb_wes200k_sig_genes_genotype_regions"

set_up_hail
set_up_pythonpath_legacy
mkdir -p ${out_dir}
python3 "${hail_script}" \
   --phenotype ${phenotype} \
   --padding 1000000 \
   --gene_table ${gene_table} \
   --final_sample_list ${final_sample_list} \
   --out_prefix ${out_prefix} 

print_update "Hail finished writing."


