#!/usr/bin/env bash
#
#$ -N permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/permute.log
#$ -e logs/permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 12
#$ -tc 10
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly bash_script="scripts/permute/_permute.sh"

readonly chr="${SGE_TASK_ID}"

# setup directories
readonly in_dir="data/permute/genes/chr${chr}"
readonly out_dir="data/permute/permutations/chr${chr}/GENE"
readonly pheno_dir="data/phenotypes"
readonly cond_dir="data/conditional/common/marker_mt"

# setup input and output paths
readonly annotation="pLoF_damaging_missense"
readonly input_path="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_chr${chr}_GENE.tsv.gz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chr${chr}_GENE"
readonly assoc_format="ukb_eur_wes_200k_maf0to5e-2_PHENO_ANNO"
readonly cond_markers="${cond_dir}/conditional_markers_chrundefined_markers.txt.gz"

# parameters for master script
readonly min_mac=4
readonly n_replicates=1000
readonly n_start_shuffle=1000
readonly n_cutoff_shuffle=10000 #10000000
readonly n_slots_saige=1
readonly n_slots_permute=1
readonly queue_saige="short.qf"
readonly queue_permute="short.qe"
readonly queue_merge="test.qc"
readonly queue_master="test.qc"
readonly n_concurrent_jobs=30
readonly iteration=1
readonly permutation_supply=0
readonly initial_top_p=10
readonly use_prs=1

# get path to true P-value and t-stats
readonly genes_path="data/permute/overview/min_mac${min_mac}/overview_genes.tsv.gz"
readonly true_p_path="data/permute/overview/min_mac${min_mac}/overview_true_p.tsv.gz"

# count how many genes to submit for the given chromosome
readonly n_genes="$( zcat ${genes_path} | grep "chr${chr}" | wc -l)"
#readonly sge_tasks="1-${n_genes}"
readonly sge_tasks="1-100"

set -x
qsub -N "_chr${chr}_permute" \
    -q "${queue_master}" \
    -pe shmem 1 \
    -t ${sge_tasks} \
    -tc ${n_concurrent_jobs} \
    "${bash_script}" \
    "${chr}" \
    "${input_path}" \
    "${out_prefix}" \
    "${pheno_dir}" \
    "${genes_path}" \
    "${true_p_path}" \
    "${min_mac}" \
    "${n_replicates}" \
    "${n_start_shuffle}" \
    "${n_cutoff_shuffle}" \
    "${n_slots_saige}" \
    "${n_slots_permute}" \
    "${queue_saige}" \
    "${queue_permute}" \
    "${queue_merge}" \
    "${queue_master}" \
    "${annotation}" \
    "${assoc_format}" \
    "${use_prs}" \
    "${cond_markers}" \
    "${iteration}" \
    "${permutation_supply}" \
    "${initial_top_p}"
set +x


