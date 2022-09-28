#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=make_genes
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/make_genes.log
#SBATCH --error=logs/make_genes.errors.log
#SBATCH --partition=test
#SBATCH --cpus-per-task 1
#SBATCH --array=21
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly curwd="${pwd}"
readonly pheno_dir="data/phenotypes"
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/permute/_make_genes.sh"

readonly chr="${SLURM_ARRAY_TASK_ID}"
readonly in_dir="data/permute/counts"
readonly out_dir="data/permute/genes/phased_only/chr${chr}"

readonly input_path="${in_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_phased_counts_chr${chr}.mt"
readonly input_type='mt'

readonly maf="maf0to5e-2"
readonly out_prefix="${out_dir}/ukb_eur_wes_200k_pLoF_damaging_missense_chr${chr}"
readonly out_type="mt"

readonly min_mac=4
readonly genes="data/permute/overview/min_mac${min_mac}/phased_only/main_genes.tsv.gz"

mkdir -p ${out_dir}

if [ -f "${genes}" ]; then
  readonly n_tasks="$( zcat ${genes} | grep -w "chr${chr}" | grep "ENSG" | wc -l)"
  readonly slurm_tasks="1-${n_tasks}"
  readonly slurm_jname="_make_genes"
  readonly slurm_project="lindgren.prj"
  readonly slurm_queue="short"
  readonly slurm_shmem="1"
  readonly jid=$( sbatch \
    --account="${project}" \
    --job-name="${jname}" \
    --output="${slum_jname}.log" \
    --error="${slurm_jname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_shmem}" \
    --array=${slurm_tasks} \
    --parsable \
    ${bash_script} \
    ${chr} \
    ${input_path} \
    ${input_type} \
    ${out_prefix} \
    ${out_type} \
    ${genes} )

else
  >&2 echo "File ${overview} does not exists! Exiting.."
fi

