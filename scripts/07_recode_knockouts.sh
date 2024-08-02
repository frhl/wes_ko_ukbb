#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=recode_knockouts
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/recode_knockouts.log
#SBATCH --error=logs/recode_knockouts.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=22

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly bash_script="scripts/_recode_knockouts.sh"
readonly merge_script="scripts/_recode_knockouts_merge.sh"
readonly rscript="scripts/_recode_knockouts.R"
readonly curwd=$(pwd)

readonly cluster=$( get_current_cluster )
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly ac_dir="data/vep/counts"
readonly ac_path="${ac_dir}/UKB.exome_array.variants.vep95.worst_csq_by_gene_canonical.original.counts.txt.gz"

readonly vep_dir="data/mt/prefilter/pp90"
readonly vep_path="${vep_dir}/ukb_wes_union_calls_200k_chr${chr}.loftee.worst_csq_by_gene_canonical.pp90.maf0_005.csqs.txt.gz"

readonly in_dir="data/knockouts/alt/pp90/encode_vcf_parallel"
#readonly out_dir="data/knockouts/alt/pp90/recoded_vep_ac/pLoF"
#readonly out_dir="data/knockouts/alt/pp90/recoded_vep_ac/damaging_missense"
readonly out_dir="data/knockouts/alt/pp90/recoded_vep_ac/pLoF_damaging_missense"
#readonly out_dir="data/knockouts/alt/pp90/recoded/test_test_test"

mkdir -p ${out_dir}

readonly chunk_size=10000

submit_recode_job()
{
  local input_path=${1}
  local out_prefix=${2}
  local jname="_c${chr}_recode"
  local lname="logs/_recode_knockouts"
  local project="lindgren.prj"
  local queue="short"
  local nslots="1"
  # get number of arrays to send
  local num_lines=$( zcat ${input_path} | wc -l )
  local num_chunks=$(( (${num_lines} + ${chunk_size} - 1) / ${chunk_size} ))
  local array_chunks="1-${num_chunks}" 
  local slurm_jid=$(sbatch \
    --account="${project}" \
    --job-name="${jname}" \
    --output="${lname}.log" \
    --error="${lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${queue}" \
    --cpus-per-task="${nslots}" \
    --array=${array_chunks} \
    --parsable \
    "${bash_script}" \
    "${rscript}" \
    "${input_path}" \
    "${vep_path}" \
    "${ac_path}" \
    "${chunk_size}" \
    "${num_chunks}" \
    "chr${chr}" \
    "${out_prefix}")
  # submit merge of chunks
  local outfile="${out_prefix}.txt.gz" 
  submit_merge_job ${out_prefix} ${outfile} ${slurm_jid}
}

submit_merge_job()
{
  echo "Submitting merge job."
  local regex_prefix="${1}"
  local outfile="${2}"
  local dependency="${3}"
  local jname="_c${chr}_recode_knockout_merge"
  local lname="logs/_recode_knockout_merge"
  local nslots="2"
  sbatch \
    --account="${project}" \
    --job-name="${jname}" \
    --output="${lname}.log" \
    --error="${lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${queue}" \
    --cpus-per-task="${nslots}" \
    --dependency="${dependency}" \
    "${merge_script}" \
    "${regex_prefix}" \
    "${num_chunks}" \
    "${outfile}"
}



#the_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_all.tsv.gz"
#the_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded.pLoF"
#submit_recode_job ${the_path} ${the_prefix}

#the_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_damaging_missense_all.tsv.gz"
#the_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded.damaging_missense"
#submit_recode_job ${the_path} ${the_prefix}

the_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_all.tsv.gz"
the_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded.pLoF_damaging_missense"
submit_recode_job ${the_path} ${the_prefix}

#the_path="${in_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense_all.tsv.gz"
#the_prefix="${out_dir}/ukb_eur_wes_200k_chr${chr}.pp90.recoded.pLoF_damaging_missense"
#submit_recode_job ${the_path} ${the_prefix}






 
