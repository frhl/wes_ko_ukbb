#!/usr/bin/env bash
#
# @description Annotate main MatrixTables with VEP results
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=debug_slurm
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/debug_slurm.log
#SBATCH --error=logs/debug_slurm.errors.log
#SBATCH --partition=epyc
#SBATCH --cpus-per-task 2
#SBATCH --array=21
#
#$ -N debug_slurm
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/debug_slurm.log
#$ -e logs/debug_slurm.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 2
#$ -q test.qc
#$ -t 21
#$ -V

source utils/bash_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="scripts/debug_slurm.py"

echo " hail mem: $(get_hail_memory)"

SECONDS=0
set_up_hail 0.2.97
set_up_pythonpath_legacy 
python3 ${hail_script} \
  --out_prefix "test"


