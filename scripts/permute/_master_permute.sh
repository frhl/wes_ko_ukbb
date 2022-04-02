#!/usr/bin/env bash
#
#
#$ -N _master_permute
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_master_permute.log
#$ -e logs/_master_permute.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q lindgren.qe
#$ -V

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh

readonly bash_script="scripts/permute/_gene_permute.sh"
readonly spa_script="scripts/permute/_gene_spa.sh"
readonly r_spa_script="scripts/permute/_get_spa_p.R"
readonly r_p_script="scripts/permute/_calc_empirical_p.R"
readonly merge_script="scripts/permute/_merge_spa.sh"

readonly chr=${1?Error: Missing argX}
readonly input_path=${2?Error: Missing argX}
readonly out_prefix=${3?Error: Missing argX}
readonly pheno_dir=${4?Error: Missing argX}
readonly true_p_path=${5?Error: Missing argX}
readonly min_mac=${6?Error: Missing argX}
readonly replicates=${7?Error: Missing argX}
readonly n_start_shuffle=${8?Error: Missing argX}
readonly n_cutoff_shuffle=${9?Error: Missing argX}
readonly n_slots_saige=${10?Error: Missing argX}
readonly n_slots_permute=${11?Error: Missing argX}
readonly tick_interval=${12?Error: Missing argX}
readonly tick_timeout=${13?Error: Missing argX}
readonly queue_saige=${14?Error: Missing argX}
readonly queue_permute=${15?Error: Missing argX}
readonly annotation=${16?Error: Missing argX}
readonly static_assoc=${17?Error: Missing argX}
readonly gene=${18?Error: Missing argX}


# get array of path to phenotype names
set_arr_phenos() {
  if [ ! -z ${pheno_dir} ]; then
    local pheno_bin="${pheno_dir}/filtered_phenotypes_binary_header.tsv"
    local pheno_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
    readarray -t arr_bin < ${pheno_bin}
    readarray -t arr_cts < ${pheno_cts}
    arr_phenos=("${arr_bin[@]}" "${arr_cts[@]}")
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

are_phenos_done() {
  echo "$( echo ${phenos_done[@]} | grep -v 0 | grep 1 | wc -l )"
}

divide_ceil() {
  local divide=${1}
  local by=${2}
  echo $(( (${divide}+${by}-1)/${by} ))
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
    echo "sleep tick ${cur_ticks} - ${n_found} / ${n_expected} files found."
    if [ ${cur_ticks} -ge ${ticks_max} ]; then
      raise_error "max ticks reached!"
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
  #local max_tasks_required=$( divide_ceil ${permutations_demand} ${replicates} )
  #local n_tasks_required=$(( ${max_tasks_required} - ${permutation_supply} ))
  local n_tasks_required=$(( (${permutations_demand} / ${replicates}) - ${permutation_supply} ))
  local sge_seed=3

  if [ ${n_tasks_required} -ge 1 ]; then

    echo "Running shuffle phase ${n_tasks_required} times."
    local tasks_lower_bound=$(( ${permutation_supply} + 1 ))
    local tasks_upper_bound=$(( ${permutation_supply} + ${n_tasks_required} ))
    local tasks_permute=${tasks_lower_bound}-${tasks_upper_bound}
    permutation_supply=$(( ${permutation_supply} + ${n_tasks_required}))

    local out_permute_success="${out_prefix}_${tasks_permute}"
    local out_permute_upper_bound="${out_prefix}_${tasks_upper_bound}"

    # if the permutation already exists. Skip submission..
    if [ ! -f "${out_permute_upper_bound}.vcf.gz" ]; then

      qsub -N "c${chr}_${gene}" \
          -t ${tasks_permute} \
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

      wait_for_files ${out_permute_success} ${n_tasks_required} ${tick_interval} ${tick_timeout}
      rm -f ${out_permute_success}*.SUCCESS

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
  local phenotype=${2}
  local pheno_saige_supply=${saige_supply[${phenotype}]}
  local n_tasks_required=$(( (${saige_demand} / ${replicates}) - ${pheno_saige_supply} ))

  if [ ${n_tasks_required} -ge 1 ]; then

    echo "Running saige for ${phenotype} with ${n_tasks_required} jobs."
    local tasks_lower_bound=$(( ${pheno_saige_supply} + 1 ))
    local tasks_upper_bound=$(( ${pheno_saige_supply} + ${n_tasks_required} ))
    local tasks_spa=${tasks_lower_bound}-${tasks_upper_bound}
    saige_supply[${phenotype}]=$(( ${pheno_saige_supply} + ${n_tasks_required}))

    local vcf_gene_spa="${out_prefix}"
    local out_gene_spa="${out_prefix}_${phenotype}"
    local out_spa_success="${out_gene_spa}_${tasks_spa}"
    local out_spa_upper_bound="${out_gene_spa}_${tasks_upper_bound}"

    local spa_name="spa_${gene}_${tasks_spa}"

    # if the SPA has already been performed. Skip it.
    if [ ! -f "${out_spa_upper_bound}.txt" ]; then

      qsub -N "${spa_name}" \
              -t ${tasks_spa} \
              -q ${queue_saige} \
              -pe shmem ${n_slots_saige} \
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

      wait_for_files ${out_spa_success} ${n_tasks_required} ${tick_interval} ${tick_timeout}
      rm -f ${out_spa_success}*.SUCCESS
    else
      echo >&2 "${out_spa_upper_bound} already exists. Skipping.."
    fi
  else
    >&2 echo "needed ${saige_demand} but already have $(( ${replicates} * ${pheno_saige_supply} )). Skipping.."
  fi

}

# use a table to lookup the true P-value from the primary analysis.
lookup_true() {

  # problematic when you are using strings that are subsets of otuer strings,
  # e.g. WHR and WHRadjBMI. Seraching for the first will result in two phenotypes.

  local phenotype=${1}
  local column=${2}

  # read file and subset to current gene/pheno/annotation
  local cur_assoc=$( echo ${static_assoc} | sed -e "s/PHENO/${phenotype}/g" | sed -e "s/ANNO/${annotation}/g") 
  local readfile=$( zcat ${true_p_path} | grep ${gene} | grep ${cur_assoc} )

  # return P-value or T-stat
  if [[ "${column}" == "p" ]]; then
    local true_value=$( echo ${readfile} | cut -d" " -f4 )
    echo ${true_value}
  elif [[ "${column}" == "t" ]]; then
    local true_value=$( echo ${readfile} | cut -d" " -f5 )
    echo ${true_value}
  else
    raise_error "Column must be t (t-statistic) or p (P-value)."
  fi

}

# aggregated saige output files into a single file
aggregate_saige() {

  local shuffles=${1}
  local phenotype=${2}
  local prefix="${out_prefix}_${phenotype}"
  local out_no_gz="${prefix}_merged.txt"
  local max_tasks=$(( (${shuffles} / ${replicates})  )) 
  rm -f "${out_no_gz}.gz"
  if [ ${shuffles} -ge 100 ]; then
    if [ ${max_tasks} -ge 1 ]; then
      for id in $(seq 1 ${max_tasks}); do
        file="${prefix}_${id}.txt"
         if [ -f ${file} ]; then
           echo ${file}
           if [ "${id}" == "1" ]; then
              cat "${file}" | head -n 1  >> "${out_no_gz}"
           fi
           cat "${file}" | tail -n +2  >> "${out_no_gz}"
         fi
       done
       gzip "${out_no_gz}"
       echo "Aggregated ${max_tasks} files to ${out_no_gz}."
    else
      raise_error "Need at least one task to submit merge"
    fi
  else
    raise_error "Invalid amount of shuffles"
  fi

}


################
# main scripts #
################

# * note: n_shuffle must be larger or equal to n_start_shuffle

set_up_rpy
set_arr_phenos
declare -A phenos_done
declare -A saige_supply
for pheno in ${arr_phenos[@]}; do phenos_done[${pheno}]=0; done
for pheno in ${arr_phenos[@]}; do saige_supply[${pheno}]=0; done

readonly log="${out_prefix}.log"
readonly success="${out_prefix}.SUCCESS"

SECONDS=0
n_shuffle=${n_start_shuffle}
permutation_supply=0
top_p=10
#testit=("WHR" "BMI")

while [ ${n_shuffle} -le ${n_cutoff_shuffle} ]; do

    shuffle_phase ${n_shuffle}

    #for phenotype in "${testit[@]}"; do
    for phenotype in "${arr_phenos[@]}"; do

      if [ ${phenos_done[${phenotype}]} == "0" ]; then

        set_arr_saige ${phenotype}
        in_gmat=${arr_saige[2]}
        in_var=${arr_saige[3]}
        gmat_bytes=$( file_size ${in_gmat} )
        var_bytes=$( file_size ${in_var} )

        if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then

          submit_saige ${n_shuffle} ${phenotype}
          saige_merged="${out_prefix}_${phenotype}_merged.txt.gz"
          aggregate_saige ${n_shuffle} ${phenotype}

          permuted_p=$( Rscript ${r_spa_script} --input_path "${saige_merged}" --select_min_p ${top_p} )
          true_p=$( lookup_true ${phenotype} "p" )

          if [ $( zcat ${saige_merged} | wc -l ) -ge 2 ]; then

            if [ $( echo "${true_p} >= ${permuted_p}" | bc ) ]; then

              phenos_done[${phenotype}]=1
              true_t=$( lookup_true ${phenotype} "t" )
              outfile="${out_prefix}_${phenotype}_empirical_p"
              empirical_p=$( Rscript ${r_p_script} --input_path "${saige_merged}" --true_tstat ${true_t} --true_p ${true_p} --out_prefix ${outfile})
              echo "${gene}\t${phenotype}\t${n_shuffle}\t${true_p}\t${permuted_p}\t${min_mac}\n" >> ${log}

            fi

          else
            phenos_done[${phenotype}]=1
            print_update "Stopped ${phenotype} x ${gene}. No markers with min_mac=${min_mac}." "${SECONDS}" 
          fi

        else
          phenos_done[${phenotype}]=1
          print_update "Stopped ${phenotype} x ${gene}. The gmat and variance ratio files are empty." "${SECONDS}" 
        fi

      fi
    done

    if [ $( are_phenos_done ) ]; then
      break
    else  
      n_shuffle=$(( ${n_shuffle} * 10 ))
      top_p=100
    fi

done

# finish up
print_update "Done! ${n_shuffle} permutations required." ${SECONDS}
mv ${log} ${success}



