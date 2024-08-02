#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh

readonly r_calc_emp_p="scripts/permute/_calc_empirical_p.R"
readonly r_get_spa_p="scripts/permute/_get_spa_p.R"
readonly rlogic="scripts/permute/_logic.R"
readonly r_get_true_p="scripts/permute/_get_true_p.R"

readonly gene=${1?Error: Missing arg1 (gene)}
readonly phenotype=${2?Error: Missing arg2 (prefix)}
readonly annotation=${3?Error: Missing arg3 (prefix)}
readonly static_assoc=${4?Error: Missing arg4 (prefix)}
readonly use_prs=${5?Error: Missing arg5 (prefix)}
readonly prs_available=${6?Error: Missing arg6 (prefix)}
readonly saige_merged=${7?Error: Missing arg7 (prefix)}
readonly top_p=${8?Error: Missing arg8 (prefix)}
readonly n_shuffle=${9?Error: Missing arg10 (prefix)}
readonly iteration=${10?Error: Missing arg11 (prefix)}
readonly out_prefix=${11?Error: Missing arg12 (prefix)}

# what phenotypes are tested?
readonly phenos="${out_prefix}.phenos"
# what phenotypes have been tested?
readonly outfile="${out_prefix}.permuted"
# file for aggregating permuted p-values
readonly pfile="${out_prefix}_${phenotype}.pvalues"


get_last_permuted_p() {
  if [ -f ${outfile} ]; then
    local last_p=$(cat ${outfile} | awk -v g="${gene}" -v p="${phenotype}" '$1==g && $2==p' | cut -f6 | tail -n1)
    if [ ! -z ${last_p} ]; then # && [ ${obs_count} -gt 0 ]; then
      echo ${last_p}
    else
      echo "NULL"
    fi
  else
    echo "NULL"
  fi
}

readonly total_phenos="$( cat ${phenos} | sort | uniq | wc -l)"
readonly tested_phenos="$( cat ${outfile} | grep -v "gene" | grep -w "OK" | sort | uniq | wc -l)"
if [ "${total_phenos}" -ne "${tested_phenos}" ]; then
  if [ -f "${saige_merged}" ]; then
    set_up_rpy
    if [ "$(zcat ${saige_merged} | wc -l)" -gt "1" ]; then
      # strange error with argparse below so needed to use commandArgs instead.
      readonly true_p=$( Rscript ${r_get_true_p} "${saige_merged}" "p" )
      if [[ "${true_p}" != "NA" ]]; then
        readonly last_p=$(get_last_permuted_p)
        readonly permuted_p=$(Rscript ${r_get_spa_p} --input_path "${saige_merged}" --select_min_p ${top_p} )
        readonly done=$(Rscript ${rlogic} --a ${true_p} --o "ge" --b ${permuted_p})
        readonly tmp_empirical_p=$(Rscript ${r_calc_emp_p} --input_path "${saige_merged}" --marker_id "${gene}" --out_prefix "${pfile}.tmp")
        echo "[msg]: ${phenotype} ${gene}. Last P = ${last_p}"
        echo "[msg]: ${phenotype} ${gene}. Permuted P = ${permuted_p}"
        echo "[msg]: ${phenotype} ${gene}. Empirical P = ${tmp_empirical_p}"
        if [[ "${last_p}" != "${permuted_p}" ]]; then
          if [ "${done}" -eq "1" ]; then
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
          readonly empirical_p="${tmp_empirical_p}"
          readonly the_status="OK (permuted_p not changing)"
          mv "${pfile}.tmp.txt.gz" "${pfile}.txt.gz"
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
    touch ${outfile}
    if [ "$( cat ${outfile} | wc -l)" -eq "0" ]; then 
      echo -e "gene\tphenotype\tprs\tn_shuffle\ttrue_p\tpermuted_p\tempirical_p\titeration\tcomment" >> ${outfile}
    fi
    echo -e "${gene}\t${phenotype}\t${prs_available}\t${n_shuffle}\t${true_p}\t${permuted_p}\t${empirical_p}\t${iteration}\t${the_status}" >> ${outfile}
  else
    raise_error "Missing saige_merged: ${saige_merged}"
  fi
else
  echo "Phenotypes have all been adequately permuted. Skipping writing _calc_p."
fi



