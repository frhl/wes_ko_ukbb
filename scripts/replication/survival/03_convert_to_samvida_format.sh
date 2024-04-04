#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=convert_to_samvida_format
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/convert_to_samvida_format.log
#SBATCH --error=logs/convert_to_samvida_format.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/replication/survival/03_convert_to_samvida_format.R"

# Directory for output files
readonly out_dir="data/survival/samvida/240226"

# Input files directory
readonly in_dir="data/survival/from_rap/240226"
readonly in_file1="${in_dir}/UKB.carrier_matrix.qced.eur.unrel.af05.pp0.90.pLoF_damaging_missense.keep.wes200k.b1of1.txt.gz"
readonly in_file2="${in_dir}/UKB.carrier_matrix.qced.eur.unrel.af05.pp0.90.pLoF_damaging_missense.keep.wes200k.b1of1.txt.gz"
readonly in_file3="${in_dir}/UKB.carrier_matrix.qced.eur.unrel.af05.pp0.90.pLoF_damaging_missense.remove.wes200k.b1of1.txt.gz"
readonly in_file4="${in_dir}/UKB.carrier_matrix.qced.eur.unrel.af05.pp0.90.pLoF_damaging_missense.remove.wes200k.b1of1.txt.gz"

# Create output directory if it doesn't exist
mkdir -p ${out_dir}

# Function to run the R script for a given input file
run_rscript() {
  local infile=$1
  local base_name=$(basename -- "$infile")
  local name_no_ext="${base_name%.*.*}"
  local outfile="${out_dir}/${name_no_ext}."

  Rscript "${rscript}" --input "${infile}" --out_prefix "${outfile}"
}

# Run the R script for each input file
set_up_rpy
run_rscript "${in_file1}"
run_rscript "${in_file2}"
run_rscript "${in_file3}"
run_rscript "${in_file4}"



