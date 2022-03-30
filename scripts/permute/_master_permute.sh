#!/usr/bin/env bash
#
#
#$ -N _gene_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_gene_permute.log
#$ -e logs/_gene_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly bash_script="scripts/permute/_gene_permute.sh"
readonly spa_script="scripts/permute/_gene_spa.sh"
readonly merge_script="scripts/permute/_merge_spa.sh"

readonly chr=${1?Error: Missing argX}
readonly input_path=${2?Error: Missing argX}
readonly out_prefix=${3?Error: Missing argX}
readonly pheno_dir=${4?Error: Missing argX}
readonly true_p_path=${5?Error: Missing argX}
readonly n_slots=${6?Error: Missing argX}
readonly start_n_shuffle=${7?Error: Missing argX}
readonly gene=${8?Error: Missing argX}
readonly annotation="pLoF"
readonly replicates=100

readonly input_path_gene=$(echo ${input_path} | sed -e "s/GENE/${gene}/g")
readonly out_prefix_gene=$(echo ${out_prefix} | sed -e "s/GENE/${gene}/g")

# todo
# remember to set new tick values!
# 

set_arr_phenos_all() {

  readonly pheno_bin="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
  readonly pheno_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
  readarray -t bin < ${pheno_bin}
  readarray -t cts < ${pheno_cts}
  declare arr_bin=${bin}
  declare arr_cts=${cts}
  declare phenos_all=("${arr_cts[@]}" "${arr_bin[@]}")

}

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
    >&2 echo "Error: ${1} is not a valid trait"
  fi

}


# return"binary" or "cts" depending on phenotype
get_trait_from_pheno() {

  local phenotype=${1}
  if [[ " ${arr_bin[*]} " =~ " ${phenotype} " ]]; then
    echo "binary"
  elif [[ " ${arr_cts[*]} " =~ " ${phenotype} " ]]; then
    echo "cts"
  else
    >&2 echo "Error: ${1} is not in binary/cts array"
  fi
}


# sleeps untill the desired file and count has been found
wait_for_files() {
  local file=${1}
  local n_expected=${2}
  local ticks_sec=${3}
  local ticks_max=${4}
  local cur_ticks=0
  local target_dir=$( dirname ${file} )
  local target_base=$( basename ${file} )
  local n_found=0
  while [ ${n_found} != ${n_expected} ]; do
    sleep ${ticks_sec}
    local n_found=$(ls -l ${target_dir} | grep ${target_base} | grep "SUCCESS" | wc -l )
    local cur_ticks=$(( ${cur_ticks} + 1 ))
    echo "waiting for ${target_base} [tick ${cur_ticks}] found ${n_found} equal to ${n_expected}"
    if [ ${cur_ticks} -ge ${ticks_max} ]; then
      >&2 echo "Error: max ticks reached!"
      exit 1
    fi
  done
  echo 1
}


shuffle_phase() {

  # note: requires $chr $input_path_gene $gene and 
  # $n_tasks_submitted as global variables
   
  local permutations_demand=${1}
  local n_tasks_required=$(( (${permutations_demand} / ${replicates}) - ${permutation_supply} ))
  local sge_seed=3

  if [ ${n_tasks_required} -ge 1 ]; then
    
    local tasks_permute=$(( ${permutation_supply} + 1 ))-$(( ${permutation_supply} + ${n_tasks_required} ))
    permutation_supply=$(( ${permutation_supply} + ${n_tasks_required}))
    local out_prefix_success="${out_prefix_gene}_${tasks_permute}"

    qsub -N "c${chr}_${gene}" \
      -t ${tasks_permute} \
      -q "short.qe" \
      -pe shmem ${n_slots} \
      ${bash_script} \
      ${chr} \
      ${input_path_gene} \
      ${out_prefix_gene} \
      ${out_prefix_success} \
      ${sge_seed} \
      ${gene} \
      ${replicates}

    wait_for_files ${out_prefix_success} ${n_tasks_required} 10s 100
    rm -f ${out_prefix_success}*.SUCCESS
  else
    >&2 echo "needed ${permutations_demand} but already have $(( ${replicates} * ${permutation_supply} )). Skipping.."
  fi

}

submit_saige() {

  # one big problem here is that SAIGE supply changes deending on phenotype.
  # should look up array to check how many have already been run

  local saige_demand=${1}
  local phenotype=${2}

  local vcf_gene_spa="${out_prefix_gene}"
  local out_gene_spa="${out_prefix_gene}_spa"

  local n_tasks_required=$(( (${saige_demand} / ${replicates}) - ${saige_supply} ))

  if [ ${n_tasks_required} -ge 1 ]; then

    local tasks_spa=$(( ${saige_supply} + 1 ))-$(( ${saige_supply} + ${n_tasks_required} ))
    saige_supply=$(( ${saige_supply} + ${n_tasks_required}))
    local out_spa_success="${out_gene_spa}_${tasks_spa}"

    local spa_name="spa_${gene}_${tasks_spa}"
    local merge_name="mrg_${gene}_${tasks_spa}"

    qsub -N "${spa_name}" \
            -t ${tasks_spa} \
            -q "short.qc" \
            -pe shmem 1 \
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

    wait_for_files ${out_spa_success} ${n_tasks_required} 10s 100
    rm -f ${out_prefix_success}*.SUCCESS

  else
    >&2 echo "needed ${permutations_demand} but already have $(( ${replicates} * ${permutation_supply} )). Skipping.."
  fi

}


get_sorted_saige_p() {
  
  local file=${1}
  local take_top=${2}

}


lookup_true_p() {

  local phenotype=${1}
  local gene=${2}
  local annotation=${3}
  local readfile=$( zcat ${true_p_path} | grep ${phenotype} | grep ${gene} | grep ${annotation})
  local true_p=$( echo ${readfile} | cut -f3 )
  echo "Looking up true P=${true_p} for ${phenotype} x ${gene} x ${annotation}"
  echo ${true_p}

}

echo "test 1"

set_arr_phenos_all
declare phenos_done=()
n_shuffle=${start_n_shuffle}
n_shuffle_cutoff=900
ermutation_supply=0
saige_supply=0

# only do untill we reach cutoff
while [ ${n_shuffle} -le ${n_shuffle_cutoff}]; do

  # while there are still phenotypes left to test
  if [ ${#phenos_all[@]} != ${#phenos_done[@]} ]; then

    shuffle_phase ${n_shuffle}

    # iterate over all phenotypes
    for phenotype in "${phenos_all[@]}"; do

      # check if phenotype is already done
      if [[ ! " ${phenotypes_done[*]} " =~ " ${phenotype} " ]]; then

        set_arr_saige ${phenotype}
        in_gmat=${arr_saige[2]}
        in_var=${arr_saige[3]}
        submit_saige ${n_shuffle} ${phenotype}

        true_p=$( lookup_true_p ${phenotype} ${gene} "${annotation}" )
        echo "${phenotype} ${gene} ${annotation}"

      fi
    break # <------
    done
    n_shuffle=$(( ${n_shuffle} * 10 )) 
  fi
  break # <------
done

echo "done"


