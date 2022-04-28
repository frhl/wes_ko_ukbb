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
readonly logic="scripts/permute/_logic.R"
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

readonly write_dir="$( dirname ${out_prefix})"
readonly task_log="${write_dir}/${gene}.log"
readonly task_log_errors="${write_dir}/${gene}.errors.log"
readonly saige_log="${write_dir}/${gene}.saige.log"
readonly saige_log_errors="${write_dir}/${gene}.saige.errors.log"


set_up_rpy
set_arr_phenos

arr_phenos=( "CC_combined" )
declare -A phenos_done
declare -A saige_supply
for pheno in ${arr_phenos[@]}; do phenos_done[${pheno}]=0; done
for pheno in ${arr_phenos[@]}; do saige_supply[${pheno}]=0; done

readonly log="${out_prefix}.log"
readonly success="${out_prefix}.SUCCESS"

SECONDS=0
iteration=0
n_shuffle=${n_start_shuffle}
permutation_supply=0
top_p=10
#testit=("WHR" "BMI")
#testit=("WHR" "EP_combined")

echo -e "gene\tphenotype\tn_shuffle\ttrue_p\tpermuted_p\tmin_mac" > ${log}


# check current P with permuted P (check if phenotypes are done)

## submit permute P depending on permutation_demand

### submit array saige array job depending on missing phenotype (Loop over phenotypes and send a phenotype-wise array job).

# if done, calculate emprical P



iteration=$(( ${iteration} + 1 ))
shuffle_phase ${n_shuffle}

#for phenotype in "${testit[@]}"; do
for phenotype in "${arr_phenos[@]}"; do

  if [ ${phenos_done[${phenotype}]} == "0" ]; then

    set_arr_saige ${phenotype}
    in_gmat=${arr_saige[2]}
    in_var=${arr_saige[3]}

    if [ -f ${in_gmat} ] && [ -f ${in_var} ]; then

      gmat_bytes=$( file_size ${in_gmat} )
      var_bytes=$( file_size ${in_var} )

      # ensure that the saige auxillary files contain data
      if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then

        submit_saige ${n_shuffle} ${phenotype}
        saige_merged="${out_prefix}_${phenotype}_merged.txt.gz"
        aggregate_saige ${n_shuffle} ${phenotype}

        permuted_p=$( Rscript ${r_spa_script} --input_path "${saige_merged}" --select_min_p ${top_p} )
        true_p=$( lookup_true ${phenotype} "p" )
        true_t=$( lookup_true ${phenotype} "t" )

        echo "${phenotype}: true_t = ${true_t}"
        echo "${phenotype}: true_p = ${true_p}"
        echo "${phenotype}: permuted_p = ${permuted_p}"

        >&2 echo "${phenotype}: true_t = ${true_t}"
        >&2 echo "${phenotype}: true_p = ${true_p}"
        >&2 echo "${phenotype}: permuted_p = ${permuted_p}"

        # Check that the saige file produces more than a single line
        # indicating more than one pseudo marker could be tested
        if [ $( zcat ${saige_merged} | wc -l ) -ge 2 ]; then

          # floating point logic with scientific notation will be handled by R to
          # check whether the true P-value is greater than the permuted P-value.
          if [ $( Rscript ${logic} --a ${true_p} --o "ge" --b ${permuted_p}) -eq "1" ]; then
            phenos_done[${phenotype}]=1
            outfile="${out_prefix}_${phenotype}_empirical_p"
            empirical_p=$( Rscript ${r_p_script} --input_path "${saige_merged}" --true_tstat ${true_t} --true_p ${true_p} --out_prefix ${outfile})
            echo -e "${gene}\t${phenotype}\t${n_shuffle}\t${true_p}\t${permuted_p}\t${min_mac}" >> ${log}
          fi
        else
          >&2 echo "setting ${phenotype} to done. (no markers, min_mac)"
          phenos_done[${phenotype}]=1
          print_update "Stopped ${phenotype} x ${gene}. No markers with min_mac=${min_mac}." "${SECONDS}" 
        fi
      else
        >&2 echo "setting ${phenotype} to done. (gmat/variance bytes)"
        phenos_done[${phenotype}]=1
        print_update "Stopped ${phenotype} x ${gene}. The gmat and/or variance ratio files are empty." "${SECONDS}" 
      fi
    else
      >&2 echo "setting ${phenotype} to done. (gmat/variance)"
      phenos_done[${phenotype}]=1
      print_update "Stopped ${phenotype} x ${gene}. The gmat and/or variance ratio does not exist." "${SECONDS}" 
    fi
  fi
done


echo "# are phenos done: ${phenos_done[@]}.# result: $( are_phenos_done )"

if [ "$( are_phenos_done )" -eq "1" ]; then
  echo "iteration ${iteration} completed! All phenotypes have been permuted accordingly."
  break
else  
  echo "itereration ${iteration} completed. Increasing shuffle to ${n_shuffle} * 10"
  n_shuffle=$(( ${n_shuffle} * 10 ))
  top_p=100
fi




