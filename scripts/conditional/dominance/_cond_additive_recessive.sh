#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly cluster=$( get_current_cluster )

readonly bash_script="scripts/conditional/dominance/_cond_additive_recessive_gene.sh"

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly in_csi=${3?Error: Missing arg3 (in_csi)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)}
readonly in_var=${5?Error: Missing arg5 (in_var)}
readonly grm_mtx=${6?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${7?Error: Missing arg7 (grm_sam)}
readonly min_mac=${8?Error: Missing arg8 (min_mac)}
readonly out_prefix=${9?Error: Missing arg9 (out_prefix)}
readonly sig_genes=${10?Error: Missing arg10 (sig_genes)}
readonly in_cond_additive_file=${11?Error: Missing arg11 (cond_additive_file)}
readonly in_cond_common_file=${12?Error: Missing arg12 (cond_common_file)}
readonly cond_annotation=${13?Error: Missing arg13 (cond_annotation)}
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

# check that the phenotypes have significant genes to run
readonly n_sig_genes=$(zcat ${sig_genes} | grep -w ${phenotype} | grep -w "chr${chr}"  | wc -l)

# Need to change CHR input depending on current task-id
readonly gmat=$(echo ${in_gmat} | sed -e "s/CHR/${chr}/g")
readonly var=$(echo ${in_var} | sed -e "s/CHR/${chr}/g")
readonly vcf=$(echo ${in_vcf} | sed -e "s/CHR/${chr}/g")
readonly csi=$(echo ${in_csi} | sed -e "s/CHR/${chr}/g")
readonly out=$(echo ${out_prefix} | sed -e "s/CHR/${chr}/g")
readonly cond_additive_file=$(echo ${in_cond_additive_file} | sed -e "s/CHR/${chr}/g")
readonly cond_common_file=$(echo ${in_cond_common_file} | sed -e "s/CHR/${chr}/g")

# Check that SAIGE step 1 has been run
readonly var_bytes=$( file_size ${var} )
readonly gmat_bytes=$( file_size ${gmat} )

if [ "${n_sig_genes}" -gt 0 ]; then
  >&2 echo "Found ${n_sig_genes} significant genes on chr${chr} for ${phenotype}."
  readonly gene_array="1-${n_sig_genes}"
  readonly slurm_jname="_c${chr}_cond_additive_recessive_gene"
  readonly slurm_lname="logs/_cond_additive_recessive_gene"
  readonly slurm_project="lindgren.prj"
  readonly slurm_queue="short"
  readonly slurm_nslots="2"
  if [ "${cluster}" = "slurm" ]; then
    sbatch \
      --account="${slurm_project}" \
      --job-name="${slurm_jname}" \
      --output="${slurm_lname}.log" \
      --error="${slurm_lname}.errors.log" \
      --chdir="$(pwd)" \
      --partition="${slurm_queue}" \
      --cpus-per-task="${slurm_nslots}" \
      --array=${gene_array} \
      --parsable \
      "${bash_script}" \
      "${phenotype}" \
      "${vcf}" \
      "${csi}" \
      "${gmat}" \
      "${var}" \
      "${grm_mtx}" \
      "${grm_sam}" \
      "${min_mac}" \
      "${out}" \
      "${sig_genes}" \
      "${cond_additive_file}" \
      "${cond_common_file}" \
      "${cond_annotation}" \
      "${chr}"
  else
    >&2 echo "${cluster} is not a valid cluster."
  fi  
else
  >&2 echo "No significant genes on chr${chr} for ${phenotype}."
fi

