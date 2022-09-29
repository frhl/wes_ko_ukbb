#!/usr/bin/env bash
#
# @description: amalgamate variants by phase to infer knockouts by genes.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=knockouts
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/knockouts.log
#SBATCH --error=logs/knockouts.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue


set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/_knockouts.sh"

readonly in_dir="data/mt/annotated"
readonly out_dir="data/knockouts/alt/only_homs"
readonly in_prefix="${in_dir}/ukb_eur_wes_union_calls_200k_chrCHR.mt"
readonly in_type="mt"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k"
readonly out_type="vcf"

# Note: ~24 slots are needed for running chr1. 
# Note: long queue may be required for chr1.
readonly tasks="21" # 1-22
readonly queue="short"
readonly project="lindgren.prj"

# should only VCF be produced?
readonly only_vcf=""

# should singletons be removed? Set to empty for FALSE
#readonly discard_prob_dosages="Y"
readonly discard_prob_dosages=""
#readonly aggr_method="collect" # either fasts or collect

# variant and sample parameters
readonly exclude="data/genes/220310_common_plofs_to_exclude.txt"
readonly af_min=""
readonly af_max=""
readonly maf_lb="0"
readonly maf_ub="5e-2"
readonly sex="both"

mkdir -p ${out_dir}

submit_knockout_job() 
{
  # I/O
  local annotation=${1}
  local nslots=${2}
  local aggr_method=${3}
  local out_prefix_csqs="${out_prefix}_chrCHR_maf${maf_lb}to${maf_ub}_${annotation/,/_}"
  local out_checkpoint="${out_prefix_csqs}_checkpoint.mt"
  # slurm specific paramters 
  local slurm_tasks="${tasks}"
  local slurm_jname="_ko_${annotation}"
  local slurm_lname="logs/_knockouts"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local slurm_nslots="${nslots}"
  readonly jid=$( sbatch \
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
    "${in_prefix}" \
    "${in_type}" \
    "${af_min}" \
    "${af_max}" \
    "${maf_lb}" \
    "${maf_ub}" \
    "${exclude}" \
    "${sex}" \
    "${annotation}" \
    "${only_vcf}" \
    "${aggr_method}" \
    "${out_prefix_csqs}" \
    "${out_type}" )
  # clean up after checkpints when
  if [ -f "${out_checkpoint}" ]; then
    rm -rf ${out_checkpoint}
  fi
}

#submit_knockout_job "synonymous" "5" "fast"
#submit_knockout_job "other_missense" "5" "fast"
#submit_knockout_job "pLoF" "5" "fast"
#submit_knockout_job "pLoF,LC" "5" "fast"
#submit_knockout_job "pLoF,LC,damaging_missense" "5" "fast"
#submit_knockout_job "damaging_missense" "5" "fast"
#

#submit_knockout_job "pLoF,damaging_missense" "24" "collect"
submit_knockout_job "pLoF,damaging_missense" "6" "only_homs"
#submit_knockout_job "pLoF" "24" "collect"
#submit_knockout_job "pLoF,damaging_missense" "6" "fast"
#submit_knockout_job "pLoF,LC" "6" "fast"
#submit_knockout_job "synonymous" "6" "fast"

#submit_knockout_job "0" "5e-2" "" "damaging_missense"
#submit_knockout_job "0" "5e-2" "" "synonymous"
#submit_knockout_job "0" "5e-2" "" "pLoF,LC,damaging_missense"

#submit_knockout_job 0 0.05 "" "pLoF"
#submit_knockout_job 0 0.05 "" "synonymous"
#§submit_knockout_job 0 0.05 "" "ptv,LC"
