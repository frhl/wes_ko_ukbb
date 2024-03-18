#!/usr/bin/env bash
#
# @description split into chunks of samples that are then pre-phased using whatshap
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=merge_prephased
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/merge_prephased.log
#SBATCH --error=logs/merge_prephased.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 5
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

# how many samples should there be in each chunk 
readonly n_split=500

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

# parameters relating to file name
readonly queue="short"
readonly samples_per_chunk=100

# save results in final directory
readonly out_dir="data/prephased/wes_union_calls"
readonly out_prefix="${out_dir}/ukb_wes_union_calls_200k_chr${chr}"

# in directories and out paths
readonly in_dir="data/prephased/wes_union_calls/chunks"
readonly in_prefix="${in_dir}/ukb_wes_union_calls_200k_chr${chr}"
readonly input_list="${in_prefix}_spc${samples_per_chunk}_${queue}.mergelist"
readonly input_super_list="${in_prefix}_spp${n_split}_${queue}.mergelist"

# save temporary uniq files here

readonly tmp="${input_list}.tmp"
readonly tmp_super="${input_super_list}.tmp"

# keep track of merged VCF

readonly out_vcf="${out_prefix}.vcf"
readonly out_vcf_gz="${out_vcf}.gz"

# keep track of scaffold (final VCF)
readonly final_vcf="${out_prefix}_scaffold"
readonly out_type="vcf"

# remove duplicates and create temporary
cat ${input_list} | sort | uniq  > ${tmp}
readonly n=$( cat ${tmp} | wc -l)
readonly n_rounded=$( echo $n | sed 's|.*|(&+500)/1000*1000|' | bc )


# this piece of code just iterates over various subsets of
# read-backed phased genotypes for 100 samples at a time, and
# combines into a larger VCF.
if [ ! -f "${out_vcf_gz}" ]; then
  module load BCFtools/1.12-GCC-10.3.0
  #module load BCFtools/1.10.2-GCC-8.3.0
  # loop over chunks of n_split
  echo "n_split=${n_split}"
  echo "n_rounded=${n_rounded}"
  # if there are more than 500 x 1000 samples them into chunks
  if [ ${n_rounded} -ge ${n_split} ]; then
    for idx_start in $(seq 1 ${n_split} ${n_rounded}); do 
      echo "Starting loop at idx ${idx_start}."
      # create temporary mergelist
      idx_end=$( echo "${idx_start}+${n_split}-1" | bc )
      tmp_idx="${tmp}_partition${idx_start}to${idx_end}"
      sed -n "${idx_start},${idx_end} p" ${tmp} > ${tmp_idx}
      # create paths to partition files
      out_idx_vcf="${in_prefix}_partition${idx_start}_${idx_end}.vcf"
      out_idx_vcf_gz="${out_idx_vcf}.gz"
      echo "Combining partition ${idx_start} to ${idx_end} in ${out_vcf}.."
      echo "${out_idx_vcf_gz}" >> ${input_super_list}
      # combine the files
      if [ ! -f "${out_idx_vcf_gz}" ]; then
        bcftools merge -l "${tmp_idx}" -Oz -o "${out_idx_vcf_gz}"
      fi
      # index files
      echo "Indexing partition ${out_idx_vcf}.."
      if [ ! -f "${out_idx_vcf_gz}.tbi" ]; then
        make_tabix "${out_idx_vcf_gz}" "tbi"
      fi
      # remove merge indexes
      rm ${tmp_idx}
    done 
  # otherwise just merge the ones we have
  else
    mv ${input_list} ${input_super_list}
  fi

  # create final VCF file
  echo "Partitions done! Merging into final file '${out_vcf_gz}'.."
  cat ${input_super_list} | sort | uniq  > ${tmp_super}
  bcftools merge -l ${tmp_super} -Oz -o "${out_vcf_gz}"
  make_tabix "${out_vcf_gz}" "tbi"
  rm ${tmp_super}
  rm ${tmp}

  # clean up files generated during loop
  for idx_start in $(seq 1 ${n_split} ${n_rounded}); do 
    idx_end=$( echo "${idx_start}+${n_split}-1" | bc )
    tmp_idx="${tmp}_partition${idx_start}to${idx_end}"
    out_idx_vcf="${in_prefix}_partition${idx_start}_${idx_end}.vcf"
    out_idx_vcf_gz="${out_idx_vcf}.gz"
    rm "${out_idx_vcf_gz}"
    rm "${out_idx_vcf_gz}.tbi"
  done
fi




