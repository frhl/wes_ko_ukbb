#!/usr/bin/env bash
#
#
#$ -N _submit_gene_spa
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_submit_gene_spa.log
#$ -e logs/_submit_gene_spa.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

readonly bash_script="scripts/permute/_gene_spa.sh"
readonly merge_script="scripts/permute/_merge_spa.sh"

readonly chr=${1?Error: Missing arg1 (chr)}
readonly phenotype=${2?Error: Missing arg2 (phenotype)}
readonly in_vcf=${3?Error: Missing arg3 (in_vcf)}
readonly in_gmat=${4?Error: Missing arg4 (in_gmat)}
readonly in_var=${5?Error: Missing arg5 (in_var)}
readonly min_mac=${6?Error: Missing arg6 (min_mac)}
readonly annotation=${7?Error: Missing arg7 (annotation)}
readonly gene_pheno_p=${8?Error: Missing arg8 (gene_pheno_p)}
readonly overview=${9?Error: Missing arg9 (overview)}
readonly gene_spa=${10?Error: Missing arg10 (gene_spa)}
readonly p_per_job=${11?Error: Missing arg11 (p_per_job)}
readonly out_prefix=${12?Error: Missing arg12 (out_prefix)}

readonly NUM=${SGE_TASK_ID}
readonly gene="$(zcat ${overview} | grep "CH" | grep "chr${chr}" | cut -f1 | sed ${NUM}'q;d' )"

readonly in_dir="$(dirname "${in_vcf}")"
readonly vcf_gene=$(echo ${in_vcf} | sed -e "s/GENE/${gene}/g")
readonly out_gene=$(echo ${out_prefix} | sed -e "s/GENE/${gene}/g")

# How many permutations were generated - used for determining permute file suffix?
readonly max_permutations="$(zcat ${overview} | grep "CH" | grep "chr${chr}" | cut -f3 | sed ${NUM}'q;d' )"
readonly max_tasks="$(( ( ${max_permutations} / ${p_per_job} ) + 1 ))"

# count jobs needed for current gene-phenotype pair
readonly n_actual_tasks=$(( $(zcat ${gene_pheno_p} | grep "${annotation}" | grep "${phenotype}" | grep ${gene} | cut -f4) + 1 ))
readonly actual_tasks="1-${n_actual_tasks}"

# check whether current gene has been tested in this phenotype
# todo: this is actually tested above..
readonly pheno_check="$( zcat ${gene_spa} | grep ${annotation} | grep ${phenotype} | grep ${gene} | wc -l )"
readonly task_limit=500


if [ $( ls ${in_dir} | grep ${gene} | wc -l) -ge 1 ]; then
  if [[ ${n_actual_tasks} -le ${task_limit} ]]; then
    if [[ ${pheno_check} -ge 0 ]]; then
      readonly merge_name="_mrg_${gene}"
      readonly spa_name="_c${chr}_${gene}"
      set -x
      # individual SPAs
      qsub -N "${spa_name}" \
          -t ${actual_tasks} \
          -q "short.qc" \
          -pe shmem 1 \
          "${bash_script}" \
          "${chr}" \
          "${vcf_gene}" \
          "${out_gene}" \
          "${in_gmat}" \
          "${in_var}" \
          "${phenotype}" \
          "${gene}" \
          "${max_tasks}" \
          "${min_mac}"
      # merge SPAs
      qsub -N ${merge_name} \
          -q "short.qc" \
          -pe shmem 1 \
          -hold_jid ${spa_name} \
          "${merge_script}" \
          "${out_gene}" \
          "${n_actual_tasks}" \
          "${max_tasks}"
      set +x
    else
       >2& echo "${gene} is not tested with ${pheotype}. Skipping.."
    fi
  else
    >&2 echo "Error: ${gene} with ${n_actual_tasks} tasks is greater than the limit of ${task_limit}! Exiting.."
  fi
else
  >&2 echo "Error: ${vcf_gene}* prefix does not exist! Exiting.."
fi





