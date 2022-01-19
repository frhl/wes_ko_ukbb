#!/usr/bin/env bash
#
# Create sparse genetic relatedness matrix using genotyped/imputed data.
#
#$ -N merge_markers
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/merge_markers.log
#$ -e logs/merge_markers.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc
# -q long.qc@@long.hge


source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly out_dir="data/prs/kinship/markers"
readonly out_prefix="${out_dir}/ukb_imp_500k_sparse_autosomes"
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

readonly hail_script="scripts/06_create_grm.py"
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"

# Generate a sequence of chromosomes to be included
chroms=$( seq 1 22 | tr '\n' ' ' )
#chroms=21

# combine markers from UKBB imputed data
if [ $( ls -1 ${out_prefix}.{bed,bim,fam} 2> /dev/null | wc -l ) -ne 3 ]; then
  set_up_hail
  set_up_pythonpath_legacy
  mkdir -p ${out_dir}
  python3 "${hail_script}" \
   --chroms ${chroms} \
   --out_prefix ${out_prefix} \
   --final_sample_list ${final_sample_list} \
   --subset_markers_by_kinship \
  conda deactivate
  print_update "Hail finished writing."
else
  print_update "${out_prefix}.bed/bim/fam already exists. Skipping"
fi



