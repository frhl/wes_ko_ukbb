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

readonly bash_script="scripts/permute/_gene_permute.sh"
readonly spa_script="scripts/permute/_gene_spa.sh"
readonly r_spa_script="scripts/permute/_get_spa_p.R"
readonly merge_script="scripts/permute/_merge_spa.sh"

readonly chr=${1?Error: Missing argX}
readonly input_path=${2?Error: Missing argX}
readonly out_prefix=${3?Error: Missing argX}
readonly pheno_dir=${4?Error: Missing argX}
readonly true_p_path=${5?Error: Missing argX}
readonly n_slots=${6?Error: Missing argX}
readonly start_n_shuffle=${7?Error: Missing argX}
readonly gene=${8?Error: Missing argX}


## consider bringing these arguments to the input
readonly annotation="pLoF"
readonly replicates=100
readonly min_mac=1
readonly static_assoc="ukb_eur_wes_200k_maf0to5e-2_PHENO_ANNO"


readonly input_path_gene=$(echo ${input_path} | sed -e "s/GENE/${gene}/g")
readonly out_prefix_gene=$(echo ${out_prefix} | sed -e "s/GENE/${gene}/g")

# todo
# remember to set new tick values!
# 

# get array of path to phenotype names
set_arr_phenos() {
  if [ ! -z ${pheno_dir} ]; then
    local pheno_bin="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
    local pheno_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
    readarray -t arr_bin < ${pheno_bin}
    readarray -t arr_cts < ${pheno_cts}
    arr_phenos=("${arr_bin[@]}" "${arr_cts[@]}")
  else
    echo "Error: global variable 'pheno_dir' has not been defined"
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
    >&2 echo "Error: ${1} is not a valid trait"
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
    >&2 echo "Error: ${1} is not in binary/cts array"
  fi
}


# sleep until the desired file and count has been found
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
      echo 0
      break
    fi
  done
  echo 1
}


# call qsub to shuffle phase depending on how many shuffles that has already
# been completed. Will wait for the qsub script to run before proceeding.
shuffle_phase() {

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

# call qsub on saige for a desired phenotype. Wait for qsub before proceeding.
submit_saige() {

  local saige_demand=${1}
  local phenotype=${2}
  local pheno_saige_supply=${saige_supply[${phenotype}]}
  local n_tasks_required=$(( (${saige_demand} / ${replicates}) - ${pheno_saige_supply} ))

  if [ ${n_tasks_required} -ge 1 ]; then

    local tasks_spa=$(( ${pheno_saige_supply} + 1 ))-$(( ${pheno_saige_supply} + ${n_tasks_required} ))
    saige_supply[${phenotype}]=$(( ${pheno_saige_supply} + ${n_tasks_required}))

    local vcf_gene_spa="${out_prefix_gene}"
    local out_gene_spa="${out_prefix_gene}_${phenotype}"
    local out_spa_success="${out_gene_spa}_${tasks_spa}"

    local spa_name="spa_${gene}_${tasks_spa}"

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
    >&2 echo "needed ${saige_demand} but already have $(( ${replicates} * ${pheno_saige_supply} )). Skipping.."
  fi

}

# use a table to lookup the true P-value from the primary analysis.
lookup_true_p() {

  # problematic when you are using strings that are subsets of otuer strings,
  # e.g. WHR and WHRadjBMI. Seraching for the first will result in two phenotypes.

  # requires global variables $assoc and $true_p_path

  local phenotype=${1}
  local gene=${2}
  local annotation=${3}

  local cur_assoc=$( echo ${static_assoc} | sed -e "s/PHENO/${phenotype}/g" | sed -e "s/ANNO/${annotation}/g") 
  local readfile=$( zcat ${true_p_path} | grep ${gene} | grep ${cur_assoc} )
  local true_p=$( echo ${readfile} | cut -d" " -f4 )

  if [ -z ${true_p} ]; then
    >&2 echo "Error: true_p could not be determined for ${phenotype} x ${gene} x ${annotation}"
  else
    echo ${true_p}
  fi

}

# aggregated saige output files into a single file
aggregate_saige() {

  local shuffles=${1}
  local phenotype=${2}
  local in_gene_spa="${out_prefix_gene}_${phenotype}"
  local out_gene_mrg="${in_gene_spa}_merged.txt"

  local max_tasks=$(( (${shuffles} / ${replicates})  )) 
  local merge_name="mrg_${gene}_${phenotype}"

  if [ ${shuffles} -ge 100 ]; then
    if [ ${max_tasks} -ge 1 ]; then

      qsub -N ${merge_name} \
          -q "short.qc" \
          -pe shmem 1 \
          "${merge_script}" \
          "${in_gene_spa}" \
          "${out_gene_mrg}" \
          "${max_tasks}"

      set_up_rpy
      wait_for_files "${out_gene_mrg}" 1 10s 30
      rm -f "${out_gene_mrg}.SUCCESS"
    else
      >&2 echo "Error: need at least one task to submit merge. Exiting.."
    fi
  else
    >&2 echo "Error: Invalid amount of shuffles. Exiting.."
  fi

}


################
# main scripts #
################

set_arr_phenos
declare -A phenos_done
declare -A saige_supply
for pheno in ${arr_phenos[@]}; do phenos_done[${pheno}]=0; done
for pheno in ${arr_phenos[@]}; do saige_supply[${pheno}]=0; done

n_shuffle=${start_n_shuffle}
n_shuffle_cutoff=900
permutation_supply=0


while [ ${n_shuffle} -le ${n_shuffle_cutoff}]; do

    shuffle_phase ${n_shuffle}

    # iterate over all phenotypes
    for phenotype in "${phenos_all[@]}"; do

      # check if phenotype is already done
      if [ ${phenos_done[${phenotype}]} == "0" ]; then

        set_arr_saige ${phenotype}
        in_gmat=${arr_saige[2]}
        in_var=${arr_saige[3]}
        submit_saige ${n_shuffle} ${phenotype}
        aggregate_saige ${n_shuffle} ${phenotype}

        saige_merged="${out_prefix_gene}_${phenotype}_merged.txt.gz"

        permuted_p=$( Rscript ${r_spa_script} --input_path "${saige_merged}" --select_min_p 100)
        true_p=$( lookup_true_p ${phenotype} ${gene} ${annotation} )

        if [ ${true_p} > $(( ${permuted_p} * 0.95 )) ]; then
          phenos_done[${phenotype}]=1
          echo "Finished ${phenotype} x ${gene} at ${n_shuffle}. P-true = ${true_p}, P-permuted[100] = ${saige_p}"
        fi
      fi

    done
    n_shuffle=$(( ${n_shuffle} * 10 )) 
done

echo "done"


