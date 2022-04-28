#!/usr/bin/env bash
#
#$ -N master_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/master_permute.log
#$ -e logs/master_permute.errors.log
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

readonly bash_script="scripts/permute/_master_permute.sh"

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
readonly sge_tasks="1-${n_genes}"

# parameters for master script
readonly min_mac=4
readonly n_replicates=10000
readonly n_start_shuffle=10000
readonly n_cutoff_shuffle=10000000
readonly n_slots_saige=1
readonly n_slots_permute=5
readonly tick_interval=30
readonly tick_timeout=800 # 10 x 400 seconds
readonly queue_saige="short.qf"
readonly queue_permute="short.qe"
readonly queue_master="short.qf"
readonly n_concurrent_jobs=5

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
    "${tick_interval}" \
    "${tick_timeout}" \
    "${queue_saige}" \
    "${queue_permute}" \
    "${annotation}" \
    "${assoc_format}"
set +x


