#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=recode_knockouts
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/recode_knockouts.log
#SBATCH --error=logs/recode_knockouts.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1
#SBATCH --requeue

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/_recode_knockouts.sh"
readonly rscript="scripts/07_recode_knockouts.R"
readonly curwd=$(pwd)

readonly cluster=$( get_current_cluster )
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly vep_dir="data/mt/prefilter/final_90_loftee"
readonly vep_path="${vep_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.csqs.txt.gz"

readonly in_dir="data/knockouts/alt/pp90/combined"
readonly out_dir="data/knockouts/alt/pp90/recoded_test"

#readonly input_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_all.tsv.gz"
#readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded"

mkdir -p ${out_dir}


submit_recode_job()
{
  local input_path=${1}
  local out_prefix=${2}
  local jname="_c${chr}_recode"
  local lname="logs/_recode_knockouts"
  local project="lindgren.prj"
  local queue="long"
  local nslots="3"
  sbatch \
    --account="${project}" \
    --job-name="${jname}" \
    --output="${lname}.log" \
    --error="${lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${queue}" \
    --cpus-per-task="${nslots}" \
    --array=${array_idx} \
    "${bash_script}" \
    "${rscript}" \
    "${input_path}" \
    "${vep_path}" \
    "${out_prefix}"
}


the_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_all.tsv.gz"
the_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded.pLoF"
submit_recode_job ${the_path} ${the_prefix}

the_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_damaging_missense_all.tsv.gz"
the_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded.damaging_missense"
submit_recode_job ${the_path} ${the_prefix}

the_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_all.tsv.gz"
the_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded.pLoF_damaging_missense"
submit_recode_job ${the_path} ${the_prefix}






 
