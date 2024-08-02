#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/hail_utils.sh
source utils/qsub_utils.sh

readonly step2_SPAtests="utils/saige/step2_SPAtests_cond.R"
readonly rscript_additive="scripts/conditional/dominance/_spa_cond_additive.R"
readonly rscript_common="scripts/conditional/common/_spa_cond_common.R"
readonly rscript_merge="scripts/conditional/combined/_brute_force_cond.R"

readonly phenotype=${1?Error: Missing arg1 (phenotype)}
readonly vcf=${2?Error: Missing arg2 (in_vcf)}
readonly csi=${3?Error: Missing arg3 (in_csi)}
readonly gmat=${4?Error: Missing arg4 (in_gmat)}
readonly var=${5?Error: Missing arg5 (in_var)}
readonly grm_mtx=${6?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${7?Error: Missing arg7 (grm_sam)}
readonly min_mac=${8?Error: Missing arg8 (min_mac)}
readonly out_prefix_chr=${9?Error: Missing arg9 (out_prefix)}
readonly sig_genes=${10?Error: Missing arg10 (markers_rare_ac)}
readonly cond_additive_file=${11?Error: Missing arg11 (cond_rare_file)}
readonly cond_common_file=${12?Error: Missing arg12 (cond_common_file)}
readonly cond_annotation=${13?Error: Missing arg13 (cond_annotation)}
readonly chr=${14?Error: Missing arg13 (cond_annotation)}

readonly index=$( get_array_task_id )
readonly gene_id=$( zcat ${sig_genes} | grep -w ${phenotype} | grep -w "chr${chr}" | sed "${index}q;d" | cut -f3 )
readonly out_prefix="${out_prefix_chr}_${gene_id}"
readonly out="${out_prefix}"

echo "Starting testing on chr${chr} for ${gene_id}."

set_up_rpy
set -x

# create a subset of common (non-coding) makers to be used in analysis
readonly out_common_markers_file="${out_prefix/CHR/${chr}}.common.markers"
if [ -f ${cond_common_file} ]; then
  Rscript "${rscript_common}" \
    --infile "${cond_common_file}" \
    --phenotype "${phenotype}" \
    --outfile "${out_common_markers_file}" \
    --pheno_col 6
fi

# create subset of additive markers to be used in analysis
readonly out_additive_markers_file="${out_prefix/CHR/${chr}}.additive.markers"
Rscript "${rscript_additive}" \
  --gene_id "${gene_id}" \
  --phenotype "${phenotype}" \
  --path_markers "${cond_additive_file}" \
  --outfile "${out_additive_markers_file}"

# merge the two subsets
readonly out_markers_file="${out_prefix/CHR/${chr}}.final.markers"
Rscript ${rscript_merge} \
  --file_common_markers ${out_common_markers_file} \
  --file_collapsed_markers ${out_additive_markers_file} \
  --file_rare_markers "n/a" \
  --outfile ${out_markers_file}

spa_test() {
  Rscript "${step2_SPAtests}" \
     --vcfFile=${vcf} \
     --vcfFileIndex=${csi} \
     --vcfField="DS" \
     --sparseGRMFile=${grm_mtx} \
     --sparseGRMSampleIDFile=${grm_sam}  \
     --chrom="chr${chr}" \
     --minMAF=0.0000001 \
     --minMAC=${min_mac} \
     --GMMATmodelFile=${gmat} \
     --varianceRatioFile=${var} \
     --SAIGEOutputFile=${out} \
     --LOCO=FALSE \
     --condition_file "${out_markers_file}" \
     && print_update "Finished saddle-point approximation for chr${chr}" ${SECONDS} \
     || raise_error "Saddle-point approximation for chr${chr} failed"
}
if [ ! -f ${out} ]; then
  set +eu
  conda deactivate
  set_up_RSAIGE
  set -eu
  spa_test
else
  >&2 echo "${out} already exists. Skipping.."
fi


