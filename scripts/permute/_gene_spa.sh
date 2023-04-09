#!/usr/bin/env bash

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly spark_dir="data/tmp/spark"
readonly step2_SPAtests="utils/saige/step2_SPAtests_cond.R"
readonly rscript="scripts/conditional/common/_spa_cond_common.R"

readonly chr=${1?Error: Missing arg1 (chr)}
readonly in_vcf=${2?Error: Missing arg2 (in_vcf)}
readonly out_gene=${3?Error: Missing arg3 (out_gene)}
readonly out_spa_success=${4?Error: Missing arg3 (out_gene)}
readonly in_gmat=${5?Error: Missing arg4 (in_gmat)}
readonly in_var=${6?Error: Missing arg5 (in_var)}
readonly grm_mtx=${7?Error: Missing arg6 (grm_mtx)}
readonly grm_sam=${8?Error: Missing arg7 (grm_sam)}
readonly phenotype=${9?Error: Missing arg6 (phenotype)}
readonly gene=${10?Error: Missing arg7 (gene)}
readonly min_mac=${11?Error: Missing arg9 (min_mac)}
readonly cond_markers=${12?Error: Missing arg10 (cond_markers)}
readonly use_cond_common=${13?Error: Missing arg11 (1 or 0 - condition on common markers?)}

readonly var_bytes=$( file_size ${in_var} )
readonly gmat_bytes=$( file_size ${in_gmat} )

readonly id=${SLURM_ARRAY_TASK_ID}
readonly vcf="${in_vcf}_${id}.vcf.gz"
readonly csi="${in_vcf}_${id}.vcf.gz.csi"
readonly out_gene_task="${out_gene}_${id}.txt"

# get conditional markers for current phenotype
if [ -f "${cond_markers}" ]; then
  readonly cond_markers_chr=$(echo ${cond_markers} | sed -e "s/CHR/${chr}/g")
  readonly out_markers="${out_gene/CHR/${chr}}.common.markers"
  if [ -f "${cond_markers_chr}" ]; then
    if [ ! -f "${out_markers}" ]; then
      set_up_rpy
      Rscript ${rscript} \
        --infile ${cond_markers_chr} \
        --phenotype ${phenotype} \
        --pheno_col 6 \
        --outfile ${out_markers}
    fi
  fi
fi



if [ ! -f ${out_gene_task} ]; then
  echo "var_bytes=${var_bytes} at ${in_var}"
  echo "gmat_bytes=${gmat_bytes} at ${in_gmat}"
  if [ ${gmat_bytes} != 0 ] && [ ${var_bytes} != 0 ]; then 
    SECONDS=0
    module purge
    #set -x
    set_up_RSAIGE
    Rscript "${step2_SPAtests}"  \
       --vcfFile=${vcf} \
       --vcfFileIndex=${csi} \
       --vcfField="DS" \
       --sparseGRMFile=${grm_mtx} \
       --sparseGRMSampleIDFile=${grm_sam}  \
       --chrom="chr${chr}" \
       --minMAF=0.0000001 \
       --minMAC=${min_mac} \
       --GMMATmodelFile=${in_gmat} \
       --varianceRatioFile=${in_var} \
       --SAIGEOutputFile=${out_gene_task} \
       --condition_file=${out_markers} \
       --LOCO=FALSE
    #set +x
    rm -f "${out_gene_task}.index"
    gzip -f ${out_gene_task}
  else
    touch ${out_file_failure}
  fi
else
  >&2 echo "${out_gene_task} already exists. Skipping.."
fi


