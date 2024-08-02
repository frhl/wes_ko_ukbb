#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=append_missing
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/append_missing.log
#SBATCH --error=logs/append_missing.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --array=1
#SBATCH --dependency="afterok:19368637"

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly cluster=$( get_current_cluster)
readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly hail_script="scripts/conditional/common/03_append_missing.py"

readonly phased_sample_list="data/phenotypes/phased_sample_list.txt"
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly samples_extract="data/phenotypes/samples/kos_not_in_imp.txt"

readonly ko_dir="data/knockouts/alt/pp90/combined"
readonly ko_path="${ko_dir}/ukb_eur_wes_200k_chr${chr}_pLoF_damaging_missense.mt"
readonly ko_type="mt"

readonly imp_dir="data/unphased/imputed/common"
readonly imp_path="${imp_dir}/ukb_imp_200k_common_chr${chr}.mt"

readonly calls_dir="data/unphased/calls/liftover"
readonly calls_path="${calls_dir}/ukb_liftover_calls_500k_chr${chr}.mt"

readonly out_dir="data/unphased/imputed/common_append_missing"
readonly out_prefix="${out_dir}/ukb_imp_200k_common_append_missing_chr${chr}"
readonly out_type="mt"


mkdir -p ${out_dir}


#if [ ! -f "${out_prefix}.mt/_SUCCESS" ]; then
  rm -rf "${out_prefix}.mt"
  set_up_hail
  set_up_pythonpath_legacy
  python3 "${hail_script}" \
       --chrom ${chr} \
       --ko_path ${ko_path} \
       --ko_type ${ko_type} \
       --imp_path ${imp_path} \
       --calls_path ${calls_path} \
       --samples_extract ${samples_extract} \
       --out_prefix ${out_prefix} \
       --out_type ${out_type} \
       && print_update "Finished filtering imputed genotypes ${out_prefix}" ${SECONDS} \
       || raise_error "Filtering imputed genotypes for for ${out_prefix} failed!"
#fi
   #--samples_qc_path ${final_sample_list} \
       #--samples_phased_path ${phased_sample_list} \
     
