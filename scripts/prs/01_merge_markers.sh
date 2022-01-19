#!/usr/bin/env bash
#
# Create sparse genetic relatedness matrix using genotyped/imputed data.
#
#$ -N merge_markers
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/merge_markers.log
#$ -e logs/merge_markers.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q long.qc@@long.hge
$ -q short.qc



source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly out_dir="data/prs/kinship/markers"
readonly out_prefix="${out_dir}/ukb_imp_eur_500k_sparse_autosomes"

readonly hail_script="scripts/prs/01_merge_markers.py"
readonly threads=$(( ${NSLOTS}-1 ))

# Generate a sequence of chromosomes to be included
chroms=$( seq 1 22 | tr '\n' ' ' )

if [ $( ls -1 ${out_prefix}.{bed,bim,fam} 2> /dev/null | wc -l ) -ne 3 ]; then
  START=${SECONDS}
  set_up_hail
  set_up_pythonpath_legacy
  mkdir -p ${out_dir}
  python3 "${hail_script}" \
   --chroms ${chroms} \
   --out_prefix ${out_prefix} \
   --ancestry "eur" \
   --subset_markers_by_kinship
  DURATION=$(( ${SECONDS}-${START} ))
  print_update "Hail finished writing after ${DURATION}"
else
  print_update "${out_prefix}.bed/bim/fam already exists. Skipping"
fi



