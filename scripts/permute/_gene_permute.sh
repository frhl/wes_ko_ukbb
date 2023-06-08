#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/vcf_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly hail_script="scripts/permute/_gene_permute.py"
readonly rscript="scripts/permute/_gene_permute.R"

readonly chr=${1?Error: Missing arg1 (phenotype)}
readonly input_path=${2?Error: Missing arg2 (input_path)}
readonly out_prefix=${3?Error: Missing arg3 (out_prefix)}
readonly out_prefix_success=${4?Error: Missing arg4 (out_prefix_success)}
readonly seed=${5?Error: Missing arg5 (seed)}
readonly gene=${6?Error: Missing arg6 (gene)}
readonly replicates=${7?Error: Missing arg7 (replicates)}
readonly cond_genotypes=${8?Error: Missing arg8 (cond_genotypes)}
readonly use_cond_common=${9?Error: Missing arg9 (use_cond_common)}

# keep track of ID
readonly id=${SLURM_ARRAY_TASK_ID}
readonly sge_seed=$(( ${id} * ${seed}))
readonly out_prefix_id="${out_prefix}_${id}"
readonly out_file_success="${out_prefix_success}_${id}.SUCCESS"

# setup common conditionion if needed
if [ "${use_cond_common}" -eq "1" ]; then
  readonly enable_cond_pipeline="YES"
else
  readonly enable_cond_pipeline=""
fi

# what phenotypes are tested and have been tested?
readonly phenos="${out_prefix}.phenos"
readonly permuted="${out_prefix}.permuted"

# only continue if phenos remain to be tested
readonly total_phenos="$( cat ${phenos} | sort | uniq | wc -l)"
readonly tested_phenos="$( cat ${permuted} | grep -v "gene" | grep -w "OK" | sort | uniq | wc -l)"

if [ "${total_phenos}" -ne "${tested_phenos}" ]; then
  if [ -f "${input_path}" ]; then
    if [ ! -f "${out_prefix_id}.vcf.gz" ]; then
      set_up_rpy
      Rscript ${rscript} \
        --chrom "chr${chr}" \
        --input_path ${input_path} \
        --input_path_cond_genotypes ${cond_genotypes} \
        --permutations ${replicates} \
        --remove_invariant_markers \
        ${enable_cond_pipeline:+--enable_cond_pipeline} \
        --out_prefix ${out_prefix_id} \
        --include_original_knockout \
        --vcf_id ${gene} \
        --seed ${sge_seed} \
        && print_update "Finished permuting phase for chr${chr}-${gene} using seed ${sge_seed}" ${SECONDS} \
        || raise_error "Permuting phase for chr${chr} failed"
      module purge
      module load BCFtools/1.12-GCC-10.3.0
      bgzip "${out_prefix_id}.vcf"
      rm -f "${out_prefix_id}.vcf"
      make_tabix "${out_prefix_id}.vcf.gz" "csi"
    else
      >&2 echo "Error: ${out_prefix_id}.vcf.bgz already exists. Skipping.."
    fi
  else
    >&2 echo "Error: ${input_path} does not exist!"
  fi
else
  >&2 echo "Skipping: All phenotypes have been adequately permuted."
fi




