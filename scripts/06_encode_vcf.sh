#!/usr/bin/env bash
#
# @description: amalgamate variants by phase to infer knockouts by genes.
# @note: check parallel/ for exporting other_missense and synonymous
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=encode_vcf
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/encode_vcf.log
#SBATCH --error=logs/encode_vcf.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly curwd=$(pwd)
readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/_encode_vcf.sh"

readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/mt/prefilter/pp90"
#readonly out_dir="data/knockouts/alt/pp90/exclude_singletons"
readonly out_dir="data/knockouts/alt/pp90/only_singletons"
readonly in_prefix="${in_dir}/ukb_wes_union_calls_200k_chrCHR.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.mt"
readonly in_type="mt"

readonly out_prefix="${out_dir}/ukb_eur_wes_200k_chrCHR"
readonly out_type="vcf"

readonly queue="short"
readonly project="lindgren.prj"

# should only VCF be produced?
readonly only_vcf=""

mkdir -p ${out_dir}


submit_encode_job() 
{
  # I/O
  local annotation=${1}
  local nslots=${2}
  local aggr_method=${3}
  local out_prefix_csqs="${out_prefix}_${annotation/,/_}"
  local out_checkpoint="${out_prefix_csqs}_checkpoint.mt"
  
  # slurm specific paramters 
  local slurm_jname="_c${chr}_ko_${annotation}"
  local slurm_lname="logs/_encode_vcf"
  local slurm_project="${project}"
  local slurm_queue="${queue}"
  local sge_queue="short.qc"
  local slurm_nslots="${nslots}"
  if [ "${cluster}" = "slurm" ]; then
    set -x
    sbatch \
      --account="${slurm_project}" \
      --job-name="${slurm_jname}" \
      --output="${slurm_lname}.log" \
      --error="${slurm_lname}.errors.log" \
      --chdir="${curwd}" \
      --partition="${slurm_queue}" \
      --cpus-per-task="${slurm_nslots}" \
      --constraint="skl-compat" \
      --array=${task_id} \
      --parsable \
      "${bash_script}" \
      "${in_prefix}" \
      "${in_type}" \
      "${annotation}" \
      "${only_vcf}" \
      "${aggr_method}" \
      "${out_prefix_csqs}" \
      "${out_type}" 
  else
    >&2 echo "${cluster} is not a valid cluster."
  fi
}

# Note: Heterozygotes/Cis are not aggregated with "fast"


#submit_encode_job "pLoF,damaging_missense" "2" "only_chets"

submit_encode_job "pLoF,damaging_missense" "2" "fast"
#submit_encode_job "pLoF,damaging_missense" "2" "fast_012"

#submit_encode_job "pLoF" "14" "collect"

#submit_encode_job "damaging_missense" "2" "fast"
#submit_encode_job "pLoF,damaging_missense" "2" "fast_012"
#submit_encode_job "damaging_missense" "24" "collect"
#submit_encode_job "pLoF" "32" "collect"
#submit_encode_job "pLoF" "2" "fast"
#submit_encode_job "other_missense" "2" "fast"
#submit_encode_job "damaging_missense" "3" "fast"
#submit_encode_job "synonymous" "3" "fast"



