#!/usr/bin/env bash
#
# @description Create sparse genetic relatedness matrix using genotyped/imputed data. Note, that this
# is generated using Weis workflow outlined here: https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html#step-0-constructing-a-sparse-grm-using-the-ld-pruned-hard-called-genotypes
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=grm
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#SBATCH --output=logs/grm.log
#SBATCH --error=logs/grm.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/bash_utils.sh
source utils/hail_utils.sh

readonly prune_script="scripts/saige/_prune_plink.sh"
readonly rv_script="scripts/saige/_create_plink_rv.sh"
readonly mrg_script="scripts/saige/_merge_plink.sh"
readonly fit_script="scripts/saige/_fit_grm.sh"
readonly liftover_script="scripts/saige/_liftover_plink.sh"

readonly spark_dir="data/tmp/spark"
readonly out_dir="data/saige/grm/input/dnanexus"
readonly prefix="${out_dir}/ukb_eur_200k_grm"
readonly final_sample_list='/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
readonly samples_keep="${out_dir}/samples.keep"

readonly curwd=$(pwd)
readonly tasks="21"
readonly project="lindgren.prj"

readonly rare_markers_per_chrom=100
readonly low_freq_markers_per_chrom=100

# get NFE samples to keep (Note, PLINK2 requires two columns)
# use python to get samples that should be kept.

cat ${final_sample_list} | grep NFE | cut -d"," -f1 | awk '{print $0,$NF}' > ${samples_keep}

mkdir -p ${out_dir}

ld_prune_plink() {
  local out_prefix=${1}
  local slurm_tasks="${tasks}"
  local slurm_jname="_prune_plink"
  local slurm_lname="logs/_prune_plink"
  local slurm_project="${project}"
  local slurm_queue="short"
  local slurm_nslots="1"
  ld_prune_plink_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array=${slurm_tasks} \
    --parsable \
    "${prune_script}" \
    "${out_prefix}" \
    "${samples_keep}" )
  echo ${ld_prune_plink_jid}
}

create_plink_vr() {
  local out_prefix=${1}
  local slurm_tasks="${tasks}"
  local slurm_jname="_plink_vr"
  local slurm_lname="logs/_create_plink_vr"
  local slurm_project="${project}"
  local slurm_queue="short"
  local slurm_nslots="1"
  readonly create_plink_vr_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --array=${slurm_tasks} \
    --parsable \
    "${rv_script}" \
    "${out_prefix}" \
    "${samples_keep}" \
    "${rare_markers_per_chrom}" \
    "${low_freq_markers_per_chrom}" )
  echo "${create_plink_vr_jid}"
}

liftover_plink() {
  local in_prefix=${1}
  local out_prefix=${2}
  local dependency=${3}
  local in_type="plink"
  local out_type="plink"
  local slurm_tasks="${tasks}"
  local slurm_jname="_liftover_plink"
  local slurm_lname="logs/_liftover_plink"
  local slurm_project="${project}"
  local slurm_queue="short"
  local slurm_nslots="1"
  readonly liftover_plink_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --array=${slurm_tasks} \
    --dependency="afterok:${dependency}" \
    --partition="${slurm_queue}" \
    --cpus-per-task="${slurm_nslots}" \
    --parsable \
    "${liftover_script}" \
    "${in_prefix}"\
    "${out_prefix}"\
    "${in_type}"\
    "${out_type}" )
  echo "${liftover_plink_jid}"
}



merge_plink() {
  local in_prefix=${1}
  local out_prefix=${2}
  local dependency=${3}
  local in_type="plink"
  local out_type="plink"
  local slurm_jname="_merge_plink"
  local slurm_lname="logs/_merge_plink"
  local slurm_project="${project}"
  local slurm_queue="short"
  local slurm_nslots="6"
  readonly merge_plink_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --dependency="afterok:${dependency}" \
    --cpus-per-task="${slurm_nslots}" \
    --parsable \
    "${mrg_script}" \
    "${in_prefix}"\
    "${in_type}"\
    "${out_prefix}"\
    "${out_type}" )
  echo "${merge_plink_jid}"
}

fit_grm() {
  local in_file=${1}
  local out_prefix=${2}
  local dependency=${3}
  local slurm_jname="_fit_grm"
  local slurm_lname="logs/_fit_grm"
  local slurm_project="${project}"
  local slurm_queue="long"
  local slurm_nslots="10"
  readonly fit_grm_jid=$( sbatch \
    --account="${slurm_project}" \
    --job-name="${slurm_jname}" \
    --output="${slurm_lname}.log" \
    --error="${slurm_lname}.errors.log" \
    --chdir="${curwd}" \
    --partition="${slurm_queue}" \
    --dependency="afterok:${dependency}" \
    --cpus-per-task="${slurm_nslots}" \
    --parsable \
    "${fit_script}" \
    "${in_file}" \
    "${out_prefix}" )
  echo "${fit_grm_jid}"
}



## create actual GRM
readonly prefix_grm="${prefix}"
readonly prefix_grm_liftover="${prefix}_grch38"
readonly prefix_combined="${prefix}_merged"
readonly prefix_fitted="${prefix}_fitted"

# step A1: create LD-pruned PLINK files
pruned_jid=$(ld_prune_plink ${prefix_grm} )

# step A2: liftover pruned plink files
liftover_pruned_jid=$( liftover_plink ${prefix_grm} ${prefix_grm_liftover} ${pruned_jid} )

# step A3: merge LD-preuned PLINK files
merged_grm_jid=$( merge_plink ${prefix_grm_liftover} ${prefix_combined} ${liftover_pruned_jid} )

# step A4: create sparse GRM
fit_grm ${prefix_combined} ${prefix_fitted} ${merged_grm_jid}


## randomly sample rare markers
readonly prefix_rv="${prefix}_rv"
readonly prefix_rv_liftover="${prefix}_grch38_rv"
readonly prefix_rv_merged="${prefix}_grch38_rv_merged"

# step B1: create rare variant plink files
rv_jid=$(create_plink_vr ${prefix_rv} )

# step B2: liftover rare variant plink files
liftover_rv_jid=$( liftover_plink ${prefix_rv} ${prefix_rv_liftover} ${rv_jid} )

# step B3: combine rare variant plink files
merge_plink ${prefix_rv_liftover} ${prefix_rv_merged} ${liftover_rv_jid}









