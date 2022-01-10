#!/usr/bin/env bash
#
# extract genotypes in regions near genes that are significant in primary analysis
#
#$ -N submit_filter_genotypes
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/submit_filter_genotypes.log
#$ -e logs/submit_filter_genotypes.errors.log
#$ -P lindgren.prjc
#$ -q test.qc
#$ -t 1-44 
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly pheno_dir="data/phenotypes"
readonly gene_dir="data/conditional/common/extract_intervals"
readonly out_dir="data/conditional/common/filter_genotypes"

readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly gene_table="${gene_dir}/211111_wes200k_saige_cts_wes_saige_sig_genes_intervals.txt"
readonly pheno_list="${pheno_dir}/UKBB_WES200k_cts_phenotypes_header.txt"
readonly filter_script="scripts/conditional/_filter_genotypes.sh"

readonly index=${SGE_TASK_ID}
readonly phenotype=$( cut -f${index} ${pheno_list} )

readonly padding=1000000
readonly min_maf=0.01
readonly min_info=0.8

mkdir -p ${out_dir}

submit_filter_genotypes_job() {
  local arg_annotation="${1}"
  local arg_out_prefix="${2}"
  set -x
  qsub -N "_filter_genotypes_${1}" \
    -q "short.qc@@short.hge" \
    -t "${SGE_TASK_ID}" \
    -pe shmem 4 \
    "${filter_script}" \
    "${phenotype}" \
    "${arg_annotation}" \
    "${gene_table}" \
    "${final_sample_list}" \
    "${arg_out_prefix}" \
    "${padding}" \
    "${min_maf}" \
    "${min_info}" \
  set +x
}

#annotation='ptv'
#out_prefix="${out_dir}/211111_intervals_${annotation}_${phenotype}"
#submit_filter_genotypes_job ${annotation} ${out_prefix}

annotation='ptv_damaging_missense'
out_prefix="${out_dir}/211111_intervals_${annotation}_${phenotype}"
submit_filter_genotypes_job ${annotation} ${out_prefix}

#annotation='synonymous'
#out_prefix="${out_dir}/211111_intervals_${annotation}_${phenotype}"
#submit_filter_genotypes_job ${annotation} ${out_prefix}



