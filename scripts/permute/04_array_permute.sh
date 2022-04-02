#!/usr/bin/env bash
#
#$ -N array_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/array_permute.log
#$ -e logs/array_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q test.qc
#$ -t 21
#$ -tc 1
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly bash_script="scripts/permute/_array_permute.sh"

readonly chr="${SGE_TASK_ID}"

# setup directories
readonly in_dir="data/permute/genes/chr${chr}"
readonly out_dir="data/permute/permutations/chr${chr}/GENE"
readonly pheno_dir="data/phenotypes"

# setup input and putput paths
readonly annotation="pLoF_damaging_missense"
readonly input_path="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_chr${chr}_GENE.tsv.gz"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chr${chr}_GENE"
readonly assoc_format="ukb_eur_wes_200k_maf0to5e-2_PHENO_ANNO"

# get path to true P-value and t-stats
readonly genes_path="data/permute/overview/overview_genes.tsv.gz"
readonly true_p_path="data/permute/overview/overview_true_p.tsv.gz"

# count how many genes to submit for the given chromosome
readonly n_genes="$( zcat ${genes_path} | grep "chr${chr}" | wc -l)"
#readonly sge_tasks="1-${n_genes}"
readonly sge_tasks=3

# parameters for master script
readonly min_mac=4
readonly n_replicates=1000
readonly n_start_shuffle=1000
readonly n_cutoff_shuffle=1000000
readonly n_slots_saige=1
readonly n_slots_permute=2
readonly tick_interval=10
readonly tick_timeout=200 # 10 x 100 seconds
readonly queue_saige="short.qf"
readonly queue_permute="short.qa"
readonly queue_master="short.qe"

set -x
qsub -N "_chr${chr}_permute" \
    -q "test.qc" \
    -pe shmem 1 \
    -t ${sge_tasks} \
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
    "${tick_interval}" \
    "${tick_timeout}" \
    "${queue_saige}" \
    "${queue_permute}" \
    "${queue_master}" \
    "${annotation}" \
    "${assoc_format}"
set +x


