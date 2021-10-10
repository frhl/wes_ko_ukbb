#!/usr/bin/env bash
#
# Create sparse genetic relatedness matrix using genotyped data.
#
#$ -N create_grm
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/create_grm.log
#$ -e logs/create_grm.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q long.qc@@long.hge

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

# paths for hail script
readonly out_dir="data/saige/grm/input/eur_ukbb"
readonly out_prefix="${out_dir}/ukb_eur_wes_sparse_markers_autosomes"
readonly hail_script="utils/subscripts/create_grm_input.py"
readonly spark_dir="data/tmp/spark"
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"

# Generate a sequence of chromosomes to be included
chroms=$( seq 1 22 | tr '\n' ' ' )

# combine markers from UKBB imputed data
if [ $( ls -1 ${out_prefix}.{bed,bim,fam} 2> /dev/null | wc -l ) -ne 3 ]; then
  set_up_hail
  set_up_pythonpath
  mkdir -p ${out_dir}
  python3 "${hail_script}" \
   --chroms ${chroms} \
   --out_prefix ${out_prefix} \
   --subset_markers_by_kinship \
   --subset_samples_by_wes200k \
   --subset_samples_by_ukbb_eur
  conda deactivate
else
  print_update "${out_prefix}.bed/bim/fam already exists. Skipping"
fi


# generate GRM using SAIGE
set_up_RSAIGE
print_update "Generating GRM from plink files.. "
Rscript "${createSparseGRM}" \
  --plinkFile=${out_prefix} \
  --nThreads=4 \
  --outputPrefix=${out_prefix} \
  --numRandomMarkerforSparseKin=1000 \
  --relatednessCutoff=0.125

