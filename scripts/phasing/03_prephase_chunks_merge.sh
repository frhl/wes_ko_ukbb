#!/usr/bin/env bash
#
# @description split into chunks of samples that are then pre-phased using whatshap
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prephase_chunks_merge
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/prephase_chunks_merge.log
#SBATCH --error=logs/prephase_chunks_merge.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 5
#SBATCH --array=21
#
#$ -N prephase_chunks_merge
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/prephase_chunks_merge.log
#$ -e logs/prephase_chunks_merge.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qc
#$ -t 21
#$ -V

set -o errexit
set -o nounset


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/vcf_utils.sh

module load BCFtools/1.12-GCC-10.3.0

# how many samples should there be in each chunk 
readonly n_split=500
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

# parameters relating to file name
readonly queue="short"
readonly samples_per_chunk=100

# directories and out paths
readonly out_dir="data/phased/wes_union_calls/prephased/chunks"
readonly out_prefix="${out_dir}/ukb_eur_wes_union_calls_200k_chr${chr}"
readonly input_list="${out_prefix}_spc${samples_per_chunk}_${queue}.mergelist"
readonly input_super_list="${out_prefix}_spp${n_split}_${queue}.mergelist"
# save temporary uniq files here
readonly tmp="${input_list}.tmp"
readonly tmp_super="${input_super_list}.tmp"
# keep track of final vcf
readonly out_vcf="${out_prefix}.vcf"
readonly out_vcf_gz="${out_vcf}.gz"

# remove duplicates and create temporary
cat ${input_list} | sort | uniq  > ${tmp}
readonly n=$( cat ${tmp} | wc -l)
readonly n_rounded=$( echo $n | sed 's|.*|(&+500)/1000*1000|' | bc )

if [ ! -f "${out_vcf_gz}" ]; then
  # loop over chunks of n_split
  for idx_start in $(seq 1 ${n_split} ${n_rounded}); do 
    
    # create temporary mergelist
    idx_end=$( echo "${idx_start}+${n_split}-1" | bc )
    tmp_idx="${tmp}_partition${idx_start}to${idx_end}"
    sed -n "${idx_start},${idx_end} p" ${tmp} > ${tmp_idx}
    
    # create paths to partition files
    out_idx_vcf="${out_prefix}_partition${idx_start}_${idx_end}.vcf"
    out_idx_vcf_gz="${out_idx_vcf}.gz"
    echo "Combining partition ${idx_start} to ${idx_end} in ${out_vcf}.."
    echo "${out_idx_vcf_gz}" >> ${input_super_list}
    
    # combine the files
    if [ ! -f "${out_idx_vcf_gz}" ]; then
      bcftools merge -l ${tmp_idx} -Oz -o "${out_idx_vcf}"
    fi
    
    # zip the files
    #if [ ! -f "${out_idx_vcf_gz}" ]; then
    #  bgzip "${out_idx_vcf}"
    #fi

    # index files
    if [ ! -f "${out_idx_vcf_gz}.tbi" ]; then
      make_tabix "${out_idx_vcf_gz}" "tbi"
    fi

    # remove merge indexes
    rm ${tmp_idx}
  done 
  
  # create final VCF file
  echo "Done! Merging into final file.."
  cat ${input_super_list} | sort | uniq  > ${tmp_super}
  bcftools merge -l ${tmp_super} -Oz -o "${out_vcf_gz}"
  make_tabix "${out_vcf_gz}" "tbi"
  
  rm ${tmp_super} ${tmp}
  echo "removing ${ouf_idx_vcf}"
fi


