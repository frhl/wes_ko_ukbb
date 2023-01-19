#!/usr/bin/env bash

set -o errexit
set -o nounset

module load BCFtools/1.12-GCC-10.3.0

source utils/qsub_utils.sh
source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly intervals=${1?Error: Missing arg1 (Gene table intervals)}
readonly out_prefix=${2?Error: Missing arg3 (out prefix)}
readonly pheno_file=${3?Error: Missing arg7 (pheno_file)}
readonly trait=${4?Error: Missing arg8 (trait)}
readonly phenotype=${5?Error: Missing arg9 (phenotype)}

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/common/06_calc_min_maf.py"
readonly checkpoint="${out_prefix}_checkpoint.mt"

calc_min_maf() 
{
  if [ ! -f "${out_prefix}_min_maf.tsv" ]; then
    python3 "${hail_script}" \
       --intervals ${intervals} \
       --out_prefix ${out_prefix} \
       --pheno_file ${pheno_file} \
       --phenotype ${phenotype} \
       --trait ${trait} \
       --min_maf_by_case_control \
       && print_update "Finished filtering imputed genotypes ${out_prefix}" ${SECONDS} \
       || raise_error "Filtering imputed genotypes for for ${out_prefix} failed!"
  else
    >%2 echo "${out_prefix}.vcf.bgz already exists. Skipping.."
  fi
}

# run analysis
set +u
set_up_hail
set -u
set_up_pythonpath_legacy
calc_min_maf

