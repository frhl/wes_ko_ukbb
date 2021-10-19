#!/usr/bin/env bash
#
# Create sparse genetic relatedness matrix using genotyped/imputed data.
#
#$ -N create_grm
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/create_grm.log
#$ -e logs/create_grm.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q long.qc@@long.hge

# -q short.qc

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"
readonly out_dir="data/saige/grm/input"
readonly out_prefix="${out_dir}/211020_ukb_wes_200k_sparse_autosomes"
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'

readonly hail_script="scripts/05_create_grm.py"
readonly threads=$(( ${NSLOTS}-1 ))
readonly createSparseGRM="/well/lindgren/flassen/software/dev/SAIGE/extdata/createSparseGRM.R"

# Generate a sequence of chromosomes to be included
chroms=$( seq 1 22 | tr '\n' ' ' )

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
   --subset_samples_by_wes200k \
   --subset_samples_by_genet_eur \
   --add_rare_variants
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




