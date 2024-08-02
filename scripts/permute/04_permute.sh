#!/usr/bin/env bash
#
# @description permute genetic phase to generate empirical P-values 
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=permute
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/permute.log
#SBATCH --error=logs/permute.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly bash_script="scripts/permute/_permute.sh"

readonly curwd=$(pwd)
readonly chr="${SLURM_ARRAY_TASK_ID}"

# setup directories
readonly in_dir="data/permute/overview"
readonly out_dir="data/permute/permutations_shuffle/chr${chr}/GENE"
readonly pheno_dir="data/phenotypes"
readonly cond_dir="data/conditional/common/markers/2024"
readonly grm_dir="data/saige/grm/input/dnanexus"

# setup input and output paths
readonly annotation="pLoF_damaging_missense"
readonly input_path="${in_dir}/sample_order.txt"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_permuted_chr${chr}_GENE"
readonly assoc_format="ukb_eur_wes_200k_PHENO_ANNO"
readonly grm_mtx="${grm_dir}/ukb_eur_200k_grm_fitted_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx"
readonly grm_sam="${grm_mtx}.sampleIDs.txt"
readonly plink_file="${grm_dir}/ukb_eur_200k_grm_grch38_rv_merged"

# The markers and actual genotypes/dosages for each sample respectively
readonly cond_markers="${cond_dir}/common_conditional.markers"
readonly cond_genotypes="${cond_dir}/common_conditional.tsv.gz"

# parameters for master script
readonly min_mac=4
readonly n_replicates=1000 # 1000
readonly n_start_shuffle=1000000 #1000
readonly n_cutoff_shuffle=1000000 #10000000
readonly n_slots_saige=1
readonly n_slots_permute=2
readonly queue_saige="short"
readonly queue_permute="short"
readonly queue_merge="short"
readonly queue_master="short"
readonly iteration=0
readonly permutation_supply=0
readonly initial_top_p=10000 # 100
readonly use_prs=1
readonly use_cond_common=1

# get path to true P-value and t-stats
readonly overview_dir="data/permute/overview"
readonly genes_path="${overview_dir}/genes_to_run_2cis_2chets.tsv.gz"
readonly genes_phenos_path="${overview_dir}/phenotypes_with_2cis_2chets.txt.gz"

# count how many genes to submit for the given chromosome
echo $( zcat ${genes_path} | grep -w "chr${chr}" )
readonly n_genes="$( zcat ${genes_path} | grep -w "chr${chr}" | wc -l)"
readonly slurm_tasks="1-${n_genes}"
readonly slurm_jname="_chr${chr}_permute"
readonly slurm_lname="logs/_permute"
readonly slurm_project="lindgren.prj"
readonly slurm_queue="${queue_master}"
readonly slurm_nslots="1"

if [ "${n_genes}" -gt "0" ]; then
  set -x
  sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array=${slurm_tasks} \
    --parsable \
    "${bash_script}" \
    "${chr}" \
    "${grm_mtx}" \
    "${grm_sam}" \
    "${input_path}" \
    "${out_prefix}" \
    "${pheno_dir}" \
    "${genes_path}" \
    "${genes_phenos_path}" \
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
    "${use_cond_common}" \
    "${cond_genotypes}" \
    "${iteration}" \
    "${permutation_supply}" \
    "${initial_top_p}"
  set +x
else
  >&2 echo "No genes to run for chr${chr}. Exiting.."
fi




