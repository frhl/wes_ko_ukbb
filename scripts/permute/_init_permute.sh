#!/usr/bin/env bash
#
#
#$ -N _init_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_init_permute.log
#$ -e logs/_init_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh

readonly this_script="scripts/permute/_init_permute.sh"
readonly bash_script="scripts/permute/_gene_permute.sh"
readonly spa_script="scripts/permute/_gene_spa.sh"
readonly merge_script="scripts/permute/_merge_spa.sh"
readonly calc_script="scripts/permute/_calc_p.sh"

readonly index=${SGE_TASK_ID}

# input arguments
readonly chr=${1?Error: Missing arg1}
readonly input_path_prelim=${2?Error: Missing arg2}
readonly out_prefix_prelim=${3?Error: Missing arg3}
readonly pheno_dir=${4?Error: Missing arg4}
readonly genes_path=${5?Error: Missing arg5} # <- new arg
readonly true_p_path=${6?Error: Missing arg6}
readonly min_mac=${7?Error: Missing arg7}
readonly replicates=${8?Error: Missing arg8}
readonly n_shuffle=${9?Error: Missing arg9}
readonly n_cutoff_shuffle=${10?Error: Missing arg9}
readonly n_slots_saige=${11?Error: Missing arg10}
readonly n_slots_permute=${12?Error: Missing arg11}
readonly tick_interval=${13?Error: Missing arg12}
readonly tick_timeout=${14?Error: Missing arg13}
readonly queue_saige=${15?Error: Missing arg14}
readonly queue_permute=${16?Error: Missing arg15}
readonly queue_merge=${17?Error: Missing arg15}
readonly queue_master=${18?Error: Missing arg15}
readonly annotation=${19?Error: Missing arg17}
readonly static_assoc=${20?Error: Missing arg18}
iteration=${21?Error: Missing arg18}
permutation_supply=${22?Error: Missing arg18}
top_p=${23?Error: Missing arg18}

# set final paths depending on gene
readonly gene="$(zcat ${genes_path} | grep "chr${chr}" | cut -f1 | sed ${index}'q;d' )"
readonly input_path=$(echo ${input_path_prelim} | sed -e "s/GENE/${gene}/g")
readonly out_prefix=$(echo ${out_prefix_prelim} | sed -e "s/GENE/${gene}/g")
readonly write_dir="$( dirname ${out_prefix})"
readonly tested_phenos="${out_prefix}.phenos"
readonly empirical_p="${out_prefix}_empirical_p.txt"

# qsub names
readonly name_shuffle="_shf_${gene}"
readonly name_saige="_spa_${gene}"
readonly name_merge="_mrg_${gene}"
readonly name_calc="_p_${gene}"
readonly name_main="_i${iteration}_${gene}"

# get logs
readonly log="${write_dir}/${gene}.log"
readonly log_saige="${write_dir}/saige.log"
readonly log_errors="${write_dir}/${gene}.errors.log"
readonly log_saige_errors="${write_dir}/saige.errors.log"

# check if permute gene even exists:
if [ ! -f "${input_path}" ]; then
  raise_error "${input_path} does not exist."
elif [ $( zcat ${input_path} | wc -l ) -le 1 ]; then
  raise_error "${input_path} only contains a header!"
fi



# create files and dirs
mkdir -p ${write_dir}
touch ${tested_phenos}

#if [ ! -z ${empirical_p} ]; then
#  echo -e "gene\tphenotype\tn_shuffle\ttrue_p\tpermuted_p\tempirical_p\tstatus" > ${empirical_p}
#fi


set_arr_phenos() {
  trait=${1}
  if [ ! -z ${pheno_dir} ]; then
    local pheno_bin="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
    local pheno_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
    readarray -t arr_bin < ${pheno_bin}
    readarray -t arr_cts < ${pheno_cts}
    if [[ "${trait}" == "both" ]]; then
      arr_phenos=("${arr_bin[@]}" "${arr_cts[@]}")
    elif [[ "${trait}" == "cts" ]]; then
      arr_phenos=("${arr_cts[@]}")
    fi
  else
    raise_error "global variable 'pheno_dir' has not been defined"
  fi
}

# set array of paths to saige null models
set_arr_saige() {
  local phenotype=${1}
  local trait=$( get_trait_from_pheno ${phenotype} )
  if [[ ( $trait == "cts" ) || ( $trait == "binary") ]]; then
    local step1_dir="data/saige/output/${trait}/step1"
    local step2_dir="data/saige/output/set/${trait}/step2"
    local in_gmat="${step1_dir}/ukb_wes_200k_${phenotype}.rda"
    local in_var="${step1_dir}/ukb_wes_200k_${phenotype}.varianceRatio.txt"
    read -a arr_saige <<< "${step1_dir} ${step2_dir} ${in_gmat} ${in_var}"
  else
    raise_error "${1} is not a valid trait"
  fi
}

# lookup if phenotype is "binary" or "cts"
get_trait_from_pheno() {
  local phenotype=${1}
  if [[ " ${arr_bin[*]} " =~ " ${phenotype} " ]]; then
    echo "binary"
  elif [[ " ${arr_cts[*]} " =~ " ${phenotype} " ]]; then
    echo "cts"
  else
    raise_error "${1} is not in binary/cts array"
  fi
}


# call qsub to shuffle phase depending on how many shuffles that has already
# been completed. Will wait for the qsub script to run before proceeding.
submit_shuffle_phase() {

  local permutations_demand=${1}
  local n_tasks_required=$(( (${permutations_demand} / ${replicates}) - ${permutation_supply} ))
  local sge_seed=3

  if [ ${n_tasks_required} -ge 1 ]; then

    echo "Running shuffle phase ${n_tasks_required} time(s)."
    local tasks_lower_bound=$(( ${permutation_supply} + 1 ))
    local tasks_upper_bound=$(( ${permutation_supply} + ${n_tasks_required} ))
    local tasks_permute=${tasks_lower_bound}-${tasks_upper_bound}
    permutation_supply=$(( ${permutation_supply} + ${n_tasks_required}))

    local out_permute_success="${out_prefix}_${tasks_permute}"
    local out_permute_upper_bound="${out_prefix}_${tasks_upper_bound}"

    # if the permutation already exists. Skip submission..
    if [ ! -f "${out_permute_upper_bound}.vcf.gz" ]; then

      qsub -N "${name_shuffle}" \
          -t ${tasks_permute} \
          -o ${log} \
          -e ${log_errors} \
          -q ${queue_permute} \
          -pe shmem ${n_slots_permute} \
          ${bash_script} \
          ${chr} \
          ${input_path} \
          ${out_prefix} \
          ${out_permute_success} \
          ${sge_seed} \
          ${gene} \
          ${replicates}

    else
      echo >&2 "${out_permute_upper_bound} already exists. Skipping.."
    fi
  else
    >&2 echo "needed ${permutations_demand} but already have $(( ${replicates} * ${permutation_supply} )). Skipping.."
  fi

}

# call qsub on saige for a desired phenotype. Wait for qsub before proceeding.
submit_saige() {

  local saige_demand=${1}
  local pheno_saige_supply=$(get_saige_supply)
  local n_tasks_required=$(( (${saige_demand} / ${replicates}) - ${pheno_saige_supply} ))

  if [ ${n_tasks_required} -ge 1 ]; then

    echo "Running saige for ${phenotype} with ${n_tasks_required} job(s)."
    local tasks_lower_bound=$(( ${pheno_saige_supply} + 1 ))
    local tasks_upper_bound=$(( ${pheno_saige_supply} + ${n_tasks_required} ))
    local tasks_spa=${tasks_lower_bound}-${tasks_upper_bound}
    #saige_supply[${phenotype}]=$(( ${pheno_saige_supply} + ${n_tasks_required}))

    local vcf_gene_spa="${out_prefix}"
    local out_gene_spa="${out_prefix}_${phenotype}"
    local out_spa_success="${out_gene_spa}_${tasks_spa}"
    local out_spa_upper_bound="${out_gene_spa}_${tasks_upper_bound}"

    local spa_name="spa_${gene}_${tasks_spa}"

    # if the SPA has already been performed. Skip it.
    if [ ! -f "${out_spa_upper_bound}.txt" ]; then
      qsub -N "${name_saige}" \
        -o ${log_saige} \
        -e ${log_saige_errors} \
        -t ${tasks_spa} \
        -q ${queue_saige} \
        -pe shmem ${n_slots_saige} \
        -hold_jid ${name_shuffle} \
        "${spa_script}" \
        "${chr}" \
        "${vcf_gene_spa}" \
        "${out_gene_spa}" \
        "${out_spa_success}" \
        "${in_gmat}" \
        "${in_var}" \
        "${phenotype}" \
        "${gene}" \
        "${min_mac}"
    else
      echo >&2 "${out_spa_upper_bound} already exists. Skipping.."
    fi
  else
    >&2 echo "needed ${saige_demand} but already have $(( ${replicates} * ${pheno_saige_supply} )). Skipping.."
  fi

}


submit_merge() {
  set -x
  qsub -N "${name_merge}" \
    -o ${log} \
    -e ${log_errors} \
    -q ${queue_merge} \
    -pe shmem 1 \
    -hold_jid "${name_saige}" \
    "${merge_script}" \
    "${n_shuffle}" \
    "${replicates}" \
    "${phenotype}" \
    "${out_prefix}" \
    "${path_merged}"
  set +x 
}

submit_calc_p() {
  echo "Calculating empirical P."
  set -x
  qsub -N "${name_calc}" \
    -o ${log} \
    -e ${log_errors} \
    -q ${queue_merge} \
    -pe shmem 1 \
    -hold_jid "${name_merge}" \
    "${calc_script}" \
    "${gene}" \
    "${phenotype}" \
    "${annotation}" \
    "${static_assoc}" \
    "${path_merged}.gz" \
    "${top_p}" \
    "${true_p_path}" \
    "${n_shuffle}" \
    "${out_prefix}"
  set +x
}

resubmit_loop() {
  echo "Resubmitting loop."
  set -x
  qsub -N "${name_main}" \
    -q "${queue_master}" \
    -pe shmem "1" \
    -t ${SGE_TASK_ID} \
    -hold_jid "${name_calc}" \
    "${this_script}" \
    "${chr}" \
    "${input_path}" \
    "${out_prefix}" \
    "${pheno_dir}" \
    "${genes_path}" \
    "${true_p_path}" \
    "${min_mac}" \
    "${replicates}" \
    "${new_n_shuffle}" \
    "${n_cutoff_shuffle}" \
    "${n_slots_saige}" \
    "${n_slots_permute}" \
    "${tick_interval}" \
    "${tick_timeout}" \
    "${queue_saige}" \
    "${queue_permute}" \
    "${queue_merge}" \
    "${queue_master}" \
    "${annotation}" \
    "${static_assoc}" \
    "${iteration}" \
    "${permutation_supply}" \
    "${new_top_p}"
  set +x
}

check_if_done() {
  local file="${out_prefix}_empirical_p.txt"
  local fcount="$(cat ${file} | grep "OK" | awk -v var="${phenotype}" '$2 == var' | wc -l)"
  if [ ${fcount} -ge "1" ]; then
    echo "1"
  else
    echo "0"
  fi
}


check_if_all_done() {
  local file="${out_prefix}_empirical_p.txt"
  local obs_count="$(cat ${file} | grep "OK" | cut -f2 | sort | uniq | wc -l)"
  local expt_count="$(cat ${tested_phenos} | sort | uniq | wc -l)"
  if [ ${obs_count} -ge 1 ]; then
    if [ ${expt_count} -eq ${obs_count} ]; then
      echo "1"
    else
      echo "0"
    fi
  else
    echo "0"
  fi
}


get_saige_supply() {
  echo "$( ls "${out_prefix}_${phenotype}_"[0-9]*.txt.gz | wc -l )"
}

SECONDS=0
do_extra_loop=0
iteration=$((${iteration} + 1))
set_arr_phenos "cts"
arr_phenos=( "Alanine_aminotransferase_residual" )

echo "Starting iteration ${iteration}"
if [ ${n_shuffle} -le ${n_cutoff_shuffle} ]; then
  if [ $(check_if_all_done) -eq "0" ]; then
    submit_shuffle_phase ${n_shuffle}
    for phenotype in "${arr_phenos[@]}"; do
      if [ $(check_if_done) -eq "0" ]; then
        set_arr_saige ${phenotype}
        in_gmat=${arr_saige[2]}
        in_var=${arr_saige[3]}
        if [ -f ${in_gmat} ] && [ -f ${in_var} ]; then
          gmat_bytes=$( file_size ${in_gmat} )
          var_bytes=$( file_size ${in_var} )
          if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then
            echo ${phenotype} >> ${tested_phenos}
            path_merged="${out_prefix}_${phenotype}_merged.txt"
            if [ ! -f ${path_merged} ]; then
              submit_saige ${n_shuffle}
              submit_merge 
            fi
            submit_calc_p
            do_extra_loop=1
          fi
       fi
      fi
    done
  fi
  if [ ${do_extra_loop} -eq "1" ]; then
    new_top_p=100
    new_n_shuffle=$(( ${n_shuffle} * 10 ))
    resubmit_loop
  else
    echo "Done! Finished all inputted phenotypes."
  fi
else
  echo "Reached cutoff (${n_cutoff_shuffle}. Last phenotype: ${phenotype}). Ending loop.."
fi




