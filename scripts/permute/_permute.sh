#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh

readonly this_script="scripts/permute/_permute.sh"
readonly bash_script="scripts/permute/_gene_permute.sh"
readonly spa_script="scripts/permute/_gene_spa.sh"
readonly merge_script="scripts/permute/_merge_permuted_spa.sh"
readonly calc_script="scripts/permute/_calc_p.sh"
readonly rscript_check_prs="scripts/saige/_check_prs_ok.R"

readonly curwd=$(pwd)
readonly project="lindgren.prj"
readonly index=${SLURM_ARRAY_TASK_ID}

# input arguments
readonly chr=${1?Error: Missing arg1}
readonly grm_mtx=${2?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${3?Error: Missing arg7 (grm_sam)}
readonly input_path=${4?Error: Missing arg2}
readonly out_prefix_prelim=${5?Error: Missing arg3}
readonly pheno_dir=${6?Error: Missing arg4}
readonly genes_path=${7?Error: Missing arg5} 
readonly genes_phenos_path=${8?Error: Missing arg5} 
readonly min_mac=${9?Error: Missing arg7 (Minimum number of knockouts)}
readonly replicates=${10?Error: Missing arg8 (Shuffles per permute submission)}
readonly n_shuffle=${11?Error: Missing arg9 (How many shuffles of of phase should be done) }
readonly n_cutoff_shuffle=${12?Error: Missing arg10 (maximum allowed number of shuffles)}
readonly n_slots_saige=${13?Error: Missing arg11 (compute slots flor saige script)}
readonly n_slots_permute=${14?Error: Missing arg12 (compute slots for permute script)}
readonly queue_saige=${15?Error: Missing arg13 (What queue should be used for SAIGE gwas)}
readonly queue_permute=${16?Error: Missing arg14 (What queue should be used for permute script)}
readonly queue_merge=${17?Error: Missing arg15 (What queue should be used for merge script)}
readonly queue_master=${18?Error: Missing arg16 (What queue should be user for master script)}
readonly annotation=${19?Error: Missing arg17 (Consequence annotation)}
readonly static_assoc=${20?Error: Missing arg18 (Prefix for association file)}
readonly use_prs=${21?Error: Missing arg19 (Should PRS be used)}
readonly cond_markers=${22?Error: Missing arg20 (List of markers to condition on)}
readonly use_cond_common=${23?Error: Missing arg20 (List of markers to condition on)}
readonly cond_genotypes=${24?Error: Missing arg21 (Genotypes/dosages for conditioning markers)}
iteration=${25?Error: Missing arg22 (What is the current iteration)}
permutation_supply=${26?Error: Missing arg23 (How many permutations have been accomplished so far)}
top_p=${27?Error: Missing arg24 (What index of the lowest P-value should be used to determien covergence (10 or 100)}

# set final paths depending on gene
readonly gene="$(zcat ${genes_path} | grep -w "chr${chr}" | cut -f1 | sed ${index}'q;d' )"
readonly out_prefix=$(echo ${out_prefix_prelim} | sed -e "s/GENE/${gene}/g")
readonly write_dir="$( dirname ${out_prefix})"
readonly tested_phenos="${out_prefix}.phenos"
readonly pheno_bin_gene="${out_prefix}.totest.phenos"
readonly status_phenos="${out_prefix}.permuted"
readonly file_cutoff="${out_prefix}.cutoff"

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

# check that path to sample order exists
if [ ! -f "${input_path}" ]; then
  raise_error "${input_path} does not exist."
fi

# check that phenotypes should actually be run for the given gene.
readonly n_phenos=$(zcat ${genes_phenos_path} | tail -n +2 | grep -w "${gene}" | wc -l)
if [ "${n_phenos}" -gt "0" ]; then
  # create files and dirs
  mkdir -p ${write_dir}
  touch ${tested_phenos}
  touch ${status_phenos}
fi

# comment out this line if you only want to run the relevant gene-trait combinations
zcat "${genes_phenos_path}" | tail -n +2 | grep -w "${gene}" | cut -f1 > "${pheno_bin_gene}"
#zcat "${genes_phenos_path}" | tail -n +2 | cut -f1 > "${pheno_bin_gene}"


set_arr_phenos() {
  trait=${1}
  if [ ! -z ${pheno_dir} ]; then
    local pheno_cts="${pheno_dir}/filtered_phenotypes_cts_manual.tsv"
    readarray -t arr_bin < ${pheno_bin_gene}
    readarray -t arr_cts < ${pheno_cts}
    if [[ "${trait}" == "both" ]]; then
      arr_phenos=("${arr_bin[@]}" "${arr_cts[@]}")
    elif [[ "${trait}" == "binary" ]]; then
      arr_phenos=("${arr_bin[@]}")
    elif [[ "${trait}" == "cts" ]]; then
      arr_phenos=("${arr_cts[@]}")
    else 
      >&2 echo "trait has not been defined properly."
    fi
  else
    raise_error "global variable 'pheno_dir' has not been defined"
  fi
}

# set array of paths to saige null models
set_arr_saige() {
  local in_prefix="ukb_wes_200k"
  local phenotype=${1}
  local trait=$( get_trait_from_pheno ${phenotype} )
  if [[ ( $trait == "cts" ) || ( $trait == "binary") ]]; then
    local step1_dir="data/saige/output/${trait}/step1"
    local step2_dir="data/saige/output/${trait}/step2"
    local in_gmat="${step1_dir}/${in_prefix}_${phenotype}.rda"
    local in_var="${step1_dir}/${in_prefix}_${phenotype}.varianceRatio.txt"
    if [ "${use_prs}" -eq "1" ]; then
      set_up_rpy # required for PRS
      local in_gmat_prs="${step1_dir}/${in_prefix}_${phenotype}_chr${chr}.rda"
      local in_var_prs="${step1_dir}/${in_prefix}_${phenotype}_chr${chr}.varianceRatio.txt"
      local prs_ok=$(Rscript ${rscript_check_prs} --phenotype ${phenotype})
      if [ -f "${in_gmat_prs/CHR/21}" ] & [ -f "${in_var_prs/CHR/21}" ] & [ "${prs_ok}" -eq "1" ]; then
        echo "Note: PRS files were found (${phenotype}). Using PRS NULL models as SAIGE input."
        local in_gmat=${in_gmat_prs}
        local in_var=${in_var_prs}
      else
        >&2 echo "Using without PRS."
      fi
    fi
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

  permutation_supply=$( get_permutation_supply )
  local permutations_demand=${1}
  local n_tasks_required=$(( (${permutations_demand} / ${replicates}) - ${permutation_supply} ))
  local sge_seed=3
  if [ ${n_tasks_required} -ge 1 ]; then
    echo "Running shuffle phase ${n_tasks_required} time(s)."
    local tasks_lower_bound=$(( ${permutation_supply} + 1 ))
    local tasks_upper_bound=$(( ${permutation_supply} + ${n_tasks_required} ))
    local tasks_permute=${tasks_lower_bound}-${tasks_upper_bound}
    #permutation_supply=$(( ${permutation_supply} + ${n_tasks_required}))
    local out_permute_success="${out_prefix}_${tasks_permute}"
    local out_permute_upper_bound="${out_prefix}_${tasks_upper_bound}"
    # if the permutation already exists. Skip submission..
    if [ ! -f "${out_permute_upper_bound}.vcf.gz" ]; then
      local slurm_tasks="${tasks_permute}"
      local slurm_jname="${name_shuffle}"
      local slurm_lname_e="${log_errors}"
      local slurm_lname_o="${log}"
      local slurm_project="${project}"
      local slurm_queue="${queue_permute}"
      local slurm_nslots="${n_slots_permute}"
      qshuffle_jid=$( sbatch \
        --account="${slurm_project}" \
        --job-name="${slurm_jname}" \
        --output="${slurm_lname_o}" \
        --error="${slurm_lname_e}" \
        --open-mode="append" \
        --chdir="${curwd}" \
        --partition="${slurm_queue}" \
        --cpus-per-task="${slurm_nslots}" \
        --array=${slurm_tasks} \
        --parsable \
        ${bash_script} \
        ${chr} \
        ${input_path} \
        ${out_prefix} \
        ${out_permute_success} \
        ${sge_seed} \
        ${gene} \
        ${replicates} \
        ${cond_genotypes} \
        ${use_cond_common} )
      echo "Submitting shuffle script (JID=${qshuffle_jid})"
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
    #local spa_name="spa_${gene}_${tasks_spa}_${phenotype}"
    
    # if the SPA has already been performed. Skip it.
    if [ ! -f "${out_spa_upper_bound}.txt" ]; then
      
      local slurm_tasks="${tasks_spa}"
      local slurm_jname="${name_saige_pheno}"
      local slurm_lname_e="${log_saige_errors}"
      local slurm_lname_o="${log_saige}"
      local slurm_project="${project}"
      local slurm_queue="${queue_saige}"
      local slurm_nslots="${n_slots_saige}"
      spa_jid=$( sbatch \
        ${qshuffle_jid:+--dependency="afterok:${qshuffle_jid}"} \
        --account="${slurm_project}" \
        --job-name="${slurm_jname}" \
        --output="${slurm_lname_o}" \
        --error="${slurm_lname_e}" \
        --chdir="${curwd}" \
        --partition="${slurm_queue}" \
        --cpus-per-task="${slurm_nslots}" \
        --array=${slurm_tasks} \
        --constraint="skl-compat" \
        --open-mode="append" \
        --parsable \
        "${spa_script}" \
        "${chr}" \
        "${vcf_gene_spa}" \
        "${out_gene_spa}" \
        "${out_spa_success}" \
        "${in_gmat}" \
        "${in_var}" \
        "${grm_mtx}" \
        "${grm_sam}" \
        "${phenotype}" \
        "${gene}" \
        "${min_mac}" \
        "${cond_markers}" \
        "${use_cond_common}")
      echo "Submitting SPA on shuffled VCFs (JID=${spa_jid})"
    else
      echo >&2 "${out_spa_upper_bound} already exists. Skipping.."
    fi
  else
    >&2 echo "needed ${saige_demand} but already have $(( ${replicates} * ${pheno_saige_supply} )). Skipping.."
  fi

}


submit_merge() {
  local slurm_jname="${name_merge_pheno}"
  local slurm_lname_e="${log_errors}"
  local slurm_lname_o="${log}"
  local slurm_project="${project}"
  local slurm_queue="${queue_merge}"
  local slurm_nslots="1"
  merge_jid=$( sbatch \
    ${spa_jid:+--dependency="afterok:${spa_jid}"} \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname_o}" \
    --error="${slurm_lname_e}" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --open-mode="append" \
    --parsable \
    "${merge_script}" \
    "${n_shuffle}" \
    "${replicates}" \
    "${phenotype}" \
    "${iteration}" \
    "${out_prefix}" \
    "${path_merged}" )
  echo "Submitting merge (JID=${merge_jid}).."
}

submit_calc_p() {
  # check if actual PRS is used in SAIGE
  local slurm_jname="${name_calc_pheno}"
  local slurm_lname_e="${log_errors}"
  local slurm_lname_o="${log}"
  local slurm_project="${project}"
  local slurm_queue="${queue_merge}"
  local slurm_nslots="1"
  calc_p_jid=$( sbatch \
    ${merge_jid:+--dependency="afterok:${merge_jid}"} \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname_o}" \
    --error="${slurm_lname_e}" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --constraint="skl-compat" \
    --open-mode="append" \
    --parsable \
    "${calc_script}" \
    "${gene}" \
    "${phenotype}" \
    "${annotation}" \
    "${static_assoc}" \
    "${use_prs}" \
    "${prs_available}" \
    "${path_merged}.gz" \
    "${top_p}" \
    "${n_shuffle}" \
    "${iteration}" \
    "${out_prefix}")
  echo "${calc_p_jid}"
}

resubmit_loop() {
  local slurm_tasks="${SLURM_ARRAY_TASK_ID}"
  local slurm_jname="${name_main}"
  local slurm_lname_e="${log_errors}"
  local slurm_lname_o="${log}"
  local slurm_project="${project}"
  local slurm_queue="${queue_master}"
  local slurm_nslots="1"
  #loop_jid=( 
  sbatch \
    --dependency="afterok:${wait_on_jids}" \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname_o}" \
    --error="${slurm_lname_e}" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array=${slurm_tasks} \
    --open-mode="append" \
    --parsable \
    "${this_script}" \
    "${chr}" \
    "${grm_mtx}" \
    "${grm_sam}" \
    "${input_path}" \
    "${out_prefix}" \
    "${pheno_dir}" \
    "${genes_path}" \
    "${genes_phenos_path}" \
    "${min_mac}" \
    "${replicates}" \
    "${new_n_shuffle}" \
    "${n_cutoff_shuffle}" \
    "${n_slots_saige}" \
    "${n_slots_permute}" \
    "${queue_saige}" \
    "${queue_permute}" \
    "${queue_merge}" \
    "${queue_master}" \
    "${annotation}" \
    "${static_assoc}" \
    "${use_prs}" \
    "${cond_markers}" \
    "${use_cond_common}" \
    "${cond_genotypes}" \
    "${iteration}" \
    "${permutation_supply}" \
    "${new_top_p}"

  echo "Re-submitted main script for another iteration!"
}

check_if_done() {
  local file="${out_prefix}.permuted"
  if [ -f ${file} ]; then
    local fcount="$(cat ${file} | grep -w "OK" | awk -v var="${phenotype}" '$2 == var' | wc -l)"
    if [ ${fcount} -ge "1" ]; then
      echo "1"
    else
      echo "0"
    fi
  else
    echo "0"
  fi
}


check_if_all_done() {
  local file=${1}
  local tested=${2}
  if [ -f "${file}" ]; then
    local obs_count="$(cat ${file} | grep -w "OK" | cut -f2 | sort | uniq | wc -l)"
    local expt_count="$(cat ${tested} | sort | uniq | wc -l)"
    if [ "${obs_count}" -ge "1" ]; then
      if [ "${expt_count}" -eq "${obs_count}" ]; then
        echo "1"
      else
        echo "0"
      fi
    else
      echo "0"
    fi
  else
    echo 0
  fi
}


# a function that counts the number of permuted VCFs created
get_permutation_supply() {
  local prefix_basename="$( basename "${out_prefix}_[0-9]*.vcf.gz$" )"
  echo "$( ls $write_dir | grep -E "${prefix_basename}" | wc -l )"
}

# A function that counts the number of SAIGE associations created
get_saige_supply() {
  local prefix_basename="$( basename "${out_prefix}_${phenotype}_[0-9]*.txt.gz" )"
  echo "$( ls $write_dir | grep -E "${prefix_basename}" | wc -l )"
}

SECONDS=0
do_extra_loop=0
iteration=$((${iteration} + 1))
wait_on_jids=""
set_arr_phenos "binary"
#arr_phenos=( "spiro_visual_impairment_and_blindness" "spiro_bronchiectasis" )
#arr_phenos=( "spiro_visual_impairment_and_blindness" "spiro_dermatitis" "spiro_asthma" )
#echo "${arr_phenos[*]}"

if [ ${n_shuffle} -le ${n_cutoff_shuffle} ]; then
  readonly is_all_done=$( check_if_all_done ${status_phenos} ${tested_phenos} )
  if [ "${is_all_done}" -eq "0" ]; then
    echo "Starting iteration ${iteration}."
    submit_shuffle_phase ${n_shuffle}
    for phenotype in "${arr_phenos[@]}"; do
      name_saige_pheno="_spa_${gene}_${phenotype}"
      name_merge_pheno="_mrg_${gene}_${phenotype}"
      name_calc_pheno="_p_${gene}_${phenotype}"
      name_calc_gene="_p_${gene}"
      echo "Testing phenotype ${phenotype} at iteration ${iteration}."
      echo ${phenotype} >> ${tested_phenos}
      if [ $(check_if_done) -eq "0" ]; then
        set_arr_saige ${phenotype}
        in_gmat=${arr_saige[2]}
        in_var=${arr_saige[3]}
        if [ -f ${in_gmat} ] && [ -f ${in_var} ]; then
          gmat_bytes=$( file_size ${in_gmat} )
          var_bytes=$( file_size ${in_var} )
          prs_available=$( echo ${in_gmat} | grep chr | wc -l)
          if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then
            path_merged="${out_prefix}_${phenotype}_merged.txt"
            if [ ! -f ${path_merged} ]; then
              echo "Running phenotype ${phenotype} with SAIGE at iteration ${iteration}"
              submit_saige ${n_shuffle}
              submit_merge 
            fi
            # only re-submit loop once all of the previous jobs (per phenotype) have completed
            do_extra_loop=1
            last_jid=$(submit_calc_p)
            wait_on_jids=$( echo "${last_jid}:${wait_on_jids}" | sed 's/:$//g' )
            echo "Waiting on the following jobs before re-initating loop: ${wait_on_jids}"
          fi
       fi
      fi
    done
  fi
  if [ "${do_extra_loop}" -eq "1" ]; then
    new_top_p=100
    new_n_shuffle=$(( ${n_shuffle} * 10 ))
    resubmit_loop
  else
    echo "Done! Finished all (selected) phenotypes."
  fi
else
  touch ${file_cutoff}
  echo "Reached cutoff (${n_cutoff_shuffle}). Ending loop.."
fi




