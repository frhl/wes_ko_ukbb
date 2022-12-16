#!/usr/bin/env bash

set -o errexit
set -o nounset

module load BCFtools/1.12-GCC-10.3.0

source utils/qsub_utils.sh
source utils/vcf_utils.sh
source utils/hail_utils.sh

readonly intervals=${1?Error: Missing arg1 (Gene table intervals)}
readonly final_sample_list=${2?Error: Missing arg2 (samples to include)}
readonly out_prefix=${3?Error: Missing arg3 (out prefix)}
readonly min_maf=${4?Error: Missing arg5 (Filter variants by min MAF)}
readonly min_info=${5?Error: Missing arg1 (Filter variants by min INFO)}
readonly missing=${6?Error: Missing arg1 (Filter variants by min INFO)}
readonly pheno_file=${7?Error: Missing arg7 (pheno_file)}
readonly trait=${8?Error: Missing arg8 (trait)}
readonly phenotype=${9?Error: Missing arg9 (phenotype)}

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/conditional/common/02_filter_genotypes.py"
readonly rstatus="${out_prefix}.running"
readonly checkpoint="${out_prefix}_checkpoint.mt"

filter_genotypes() 
{
  >&2 echo "Running if statement."
  if [ ! -f "${out_prefix}.vcf.bgz" ]; then
    >&2 echo "Python will be used now:"
    touch ${rstatus} 
    python3 "${hail_script}" \
       --min_maf ${min_maf} \
       --min_info ${min_info} \
       --missing ${missing} \
       --intervals ${intervals} \
       --extract ${final_sample_list} \
       --out_prefix ${out_prefix} \
       --pheno_file ${pheno_file} \
       --phenotype ${phenotype} \
       --trait ${trait} \
       --checkpoint \
       --min_maf_by_case_control \
       && print_update "Finished filtering imputed genotypes ${out_prefix}" ${SECONDS} \
       || raise_error "Filtering imputed genotypes for for ${out_prefix} failed!"
  rm ${rstatus}
  rm -f ${checkpoint}
  else
    >%2 echo "${out_prefix}.vcf.bgz already exists. Skipping.."
  fi
}

# run analysis
set_up_hail
set_up_pythonpath_legacy
>&2 echo "starting analysis for ${phenotype}.."
filter_genotypes
make_tabix "${out_prefix}.vcf.bgz" "csi"
>&2 echo "All done for ${phenotype}"

