#!/usr/bin/env bash
#
#$ -N _calc_p
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_calc_p.log
#$ -e logs/_calc_p.errors.log
#$ -V

set -o errexit
set -o nounset

source utils/bash_utils.sh


readonly r_calc_emp_p="scripts/permute/_calc_empirical_p.R"
readonly r_get_spa_p="scripts/permute/_get_spa_p.R"
readonly rlogic="scripts/permute/_logic.R"

readonly gene=${1?Error: Missing arg1 (prefix)}
readonly phenotype=${2?Error: Missing arg1 (prefix)}
readonly annotation=${3?Error: Missing arg1 (prefix)}
readonly static_assoc=${4?Error: Missing arg1 (prefix)}
readonly use_prs=${5?Error: Missing arg1 (prefix)}
readonly prs_available=${6?Error: Missing arg1 (prefix)}
readonly saige_merged=${7?Error: Missing arg1 (prefix)}
readonly top_p=${8?Error: Missing arg1 (prefix)}
readonly true_p_path=${9?Error: Missing arg1 (prefix)}
readonly n_shuffle=${10?Error: Missing arg1 (prefix)}
readonly out_prefix=${11?Error: Missing arg1 (prefix)}

readonly outfile="${out_prefix}.permuted"
readonly pfile="${out_prefix}_${phenotype}.pvalues"

# use a table to lookup the true P-value from the primary analysis.
# problematic when you are using strings that are subsets of otuer strings,
# e.g. WHR and WHRadjBMI. Seraching for the first will result in two phenotypes. 
lookup_true() {
  local column=${1}
  # read file and subset to current gene/pheno/annotation
  if [ -f ${true_p_path} ]; then
    if [ "${use_prs}" -eq "1" ] & [ "${prs_available}" -eq "1" ]; then
      local cur_assoc=$( echo ${static_assoc} | sed -e "s/PHENO/${phenotype}/g" | sed -e "s/ANNO/${annotation}/g")
      local readfile=$( zcat ${true_p_path} | grep ${gene} | grep ${cur_assoc} | grep "locoprs" )
      local lines=$( zcat ${true_p_path} | grep ${gene} | grep ${cur_assoc} | grep "locoprs" | wc -l)
    else
      local cur_assoc=$( echo ${static_assoc} | sed -e "s/PHENO/${phenotype}/g" | sed -e "s/ANNO/${annotation}/g")
      local readfile=$( zcat ${true_p_path} | grep ${gene} | grep ${cur_assoc} | grep -v "locoprs" )
      local lines=$( zcat ${true_p_path} | grep ${gene} | grep ${cur_assoc} | grep -v "locoprs" | wc -l)
    fi
     # return P-value or T-stat
    if [ "${lines}" -eq 1 ]; then
      if [[ "${column}" == "p" ]]; then
        local true_value=$( echo ${readfile} | cut -d" " -f4 )
        echo ${true_value}
      elif [[ "${column}" == "t" ]]; then
        local true_value=$( echo ${readfile} | cut -d" " -f5 )
        echo ${true_value}
      else
        raise_error "Column must be t (t-statistic) or p (P-value)."
      fi
    elif [ "${lines}" -le 1 ]; then
      echo "NA"
      >&2 echo "Lookup true P-value failed for ${gene} at ${cur_assoc} (No lines! true_p does not exist.)"
    else
      echo "NA"
      >&2 echo "Lookup true P-value failed for ${gene} at ${cur_assoc} (Too many lines! Lines=${lines}.)"
   fi
 else
   raise_error "Overview path (${true_p_path}) does not exist!"
 fi
}

get_last_permuted_p() {
  if [ -f ${outfile} ]; then
    local last_p=$(cat ${outfile} | awk -v g="${gene}" -v p="${phenotype}" '$1==g && $2==p' | cut -f5 | tail -n1)
    #local obs_count=$(cat ${outfile} | awk -v g="${gene}" -v p="${phenotype}" '$1==g && $2==p' | cut -f5 | tail -n1 | wc -l)
    if [ ! -z ${last_p} ]; then # && [ ${obs_count} -gt 0 ]; then
      echo ${last_p}
    else
      echo "NULL"
    fi
  else
    echo "NULL"
  fi
}

if [ -f ${saige_merged} ]; then
  readonly true_p=$( lookup_true "p" )
  readonly true_t=$( lookup_true "t" )
  if [ $(zcat ${saige_merged} | wc -l) -gt 1 ]; then
    if [[ "${true_p}" != "NA" ]]; then
      set_up_rpy
      readonly last_p=$(get_last_permuted_p)
      readonly permuted_p=$(Rscript ${r_get_spa_p} --input_path "${saige_merged}" --select_min_p ${top_p} )
      readonly done=$(Rscript ${rlogic} --a ${true_p} --o "ge" --b ${permuted_p})
      readonly tmp_empirical_p=$(
          Rscript ${r_calc_emp_p} \
            --input_path "${saige_merged}" \
            --true_tstat ${true_t} \
            --true_p ${true_p} \
            --out_prefix "${pfile}.tmp")
      echo "[msg]: ${phenotype} ${gene}. Last P = ${last_p}"
      echo "[msg]: ${phenotype} ${gene}. Permuted P = ${permuted_p}"
      echo "[msg]: ${phenotype} ${gene}. Empirical P = ${tmp_empirical_p}"
      if [[ "${last_p}" != "${permuted_p}" ]]; then
        if [ ${done} -eq 1 ]; then
          readonly the_status="OK"
          readonly empirical_p=${tmp_empirical_p}
          mv "${pfile}.tmp.txt.gz" "${pfile}.txt.gz"
        else
          readonly empirical_p="NA"
          readonly the_status="NA"
        fi
      else
        # P-value is not changing which
        # should result in termination of script
        readonly empirical_p="NA"
        readonly the_status="OK (permuted_p not changing)"
      fi
    else
      # gene not present in saige 
      # primary analysis for phenotype
      readonly permuted_p="NA"
      readonly empirical_p="NA"
      readonly the_status="OK (true_p invalid)"
    fi
  else
    # min mac causing saige_merged to
    # to have only a single line
    readonly permuted_p="NA"
    readonly empirical_p="NA"
    readonly the_status="OK (min_mac)"
  fi
  echo -e "${gene}\t${phenotype}\t${prs_available}\t${n_shuffle}\t${true_p}\t${permuted_p}\t${empirical_p}\t${the_status}" >> ${outfile}
else
  raise_error "Missing saige_merged: ${saige_merged}"
fi




