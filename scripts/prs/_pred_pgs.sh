#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly r_script=${1?Error: Missing arg1 (r_script)}
readonly pred=${2?Error: Missing arg3 (prediction file)}
readonly ldsc=${3?Error: Missing arg2 (ld_matrix)}
readonly ld_dir=${4?Error: Missing arg2 (ld_matrix)}
readonly path_betas=${5?Error: Missing arg2 (betas)}
readonly impute=${6?Error: Missing arg2 (impute)}
readonly prefix=${7?Error: Missing arg8 (prefix)}

readonly cluster=$( get_current_cluster )
readonly index=$( get_array_task_id )
readonly chr=$( get_chr ${index} )

readonly pred_chr=$(echo ${pred} | sed -e "s/CHR/${chr}/g")
readonly out_prefix_chr=$(echo ${prefix} | sed -e "s/CHR/${chr}/g")

readonly tmp_bfile="${out_prefix_chr}.bfile"
readonly tmp_bk="${tmp_bfile}.bk"
readonly tmp_rds="${tmp_bfile}.rds"

if [ ! -f "${out_prefix_chr}.txt.gz" ]; then
  set_up_ldpred2
  Rscript "${r_script}" \
      --chrom "chr${chr}" \
      --pred "${pred_chr}" \
      --impute "${impute}" \
      --tmp_bfile "${tmp_bfile}" \
      --out_prefix "${out_prefix_chr}" \
      --path_weights ${path_betas} 
  # always remove temporary bk files as these
  # tend to become extremely large (In the magnitude of terrabytes)  
  rm ${tmp_bk} ${tmp_rds}
else
  echo "Note: ${out_prefix_chr} already exists. Skipping.."
fi

