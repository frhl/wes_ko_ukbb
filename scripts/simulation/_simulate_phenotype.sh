#!/usr/bin/env bash
#
#$ -N _simulate_phenotype
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/_simulate_phenotype.log
#$ -e logs/_simulate_phenotype.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge
#$ -t 1
#$ -V

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/bash_utils.sh

readonly hail_script="scripts/simulation/01_simulate_phenotype.py"
readonly rscript="scripts/simulation/01_simulate_phenotype.R"
readonly spark_dir="data/tmp/spark_dir"

readonly in_prefix=${1?Error: Missing arg1 (phenotype)}
readonly in_type=${2?Error: Missing arg2 (in_vcf)}
readonly h2_co=${3?Error: Missing arg3 ()}
readonly pi_co=${4?Error: Missing arg3 ()}
readonly h2_ko=${5?Error: Missing arg3 ()}
readonly pi_ko=${6?Error: Missing arg3 ()}
readonly K=${7?Error: Missing arg3 ()}
readonly seed=${8?Error: Missing arg3 ()}
readonly out_prefix=${9?Error: Missing arg3 ()}

readonly out_sge_prefix="${out_prefix}_${SGE_TASK_ID}"
readonly sge_seed=$(( ${SGE_TASK_ID} * ${seed}))

mkdir -p ${spark_dir}

SECONDS=0
set_up_hail
set_up_pythonpath_legacy
set -x
python3 "${hail_script}" \
   --in_prefix "${in_prefix}"\
   --in_type "${in_type}" \
   --h2_co ${h2_co} \
   --pi_co ${pi_co} \
   --h2_ko ${h2_ko} \
   --pi_ko ${pi_ko} \
   --K ${K} \
   --prune_hom_alt 0.98 \
   --max_maf 0.05 \
   --seed ${sge_seed} \
   --out_prefix "${out_sge_prefix}" \
   && print_update "Finished simulating phenotypes for ${in_prefix}" ${SECONDS} \
   || raise_error "Simulating phenotypes for ${in_prefix} failed"
set +x

