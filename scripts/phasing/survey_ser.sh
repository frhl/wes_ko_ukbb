#!/usr/bin/env bash
#
#$ -N wes_union_call_ser
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/wes_union_call_ser.log
#$ -e logs/wes_union_call_ser.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hga
#$ -t 20-22
#$ -V

module load BCFtools/1.12-GCC-10.3.0

readonly rscript="scripts/phasing/ser.R"

readonly in_dir="data/phased/wes/union_call"
readonly in="${in_dir}/"
readonly out_dir="data/switch"
readonly out="${out_dir}/ukb_wes_union_calls_200k_chr${chr}.txt"

readonly chr="${SGE_TASK_ID}"

module purge
module load R
Rscript ${rscript} \
    --switch-error-file ${ser_file} \
    --variant-file ${var_file} \
    --outfile ${out_file}


