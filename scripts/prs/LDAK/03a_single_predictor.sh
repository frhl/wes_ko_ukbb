#!/usr/bin/env bash
#
#$ -N calculate_kinship
#$ -wd /well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/calculate_kinship.log
#$ -e logs/calculate_kinship.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q long.qc@@long.hga
#$ -t 21
#$ -V


readonly in_dir="data/saige/grm/input"
readonly out_dir="data/prs/kinship/matrix"

readonly ldak_dir="/well/lindgren/flassen/software/ldak"
readonly ldak="${ldak_dir}/ldak5.2.linux"

readonly in_file="${in_dir}/211102_long_ukb_wes_200k_sparse_autosomes"
readonly out_thin="${out_dir}/thin"
readonly out_kins="${out_dir}/kins"

mkdir -p ${out_dir}


# thinning predictors
if [ ! -f "${out_thin}.weights" ]; then
  set -x
  ${ldak} \
    --thin "${out_thin}" \
    --bfile "${in_file}" \
    --window-prune .98 \
    --window-kb 100
  set +x
  awk < "${out_thin}.in" '{print $1, 1}' > "${out_thin}.weights"
else
  echo "${out_thin}.weights file found. Skipping thinning."
fi

# calc kinship matrix assuming LDAK-thin model
if [ ! -f ${out_kins} ]; then
  set -x 
  ${ldak} \
    --calc-kins-direct "${out_kins}" \
    --bfile "${in_file}" \
    --weights "${out_thin}.weights" \
    --power -0.25
  set +x
fi

#bash ${ldak} \
#  --cut-kins "${out}" \
#  --partition-length "500000" \
#  --bfile "${in_file}"

#bash "${ldak}" \
#  --calc-kins ${out} \
#  --partition 1 \
#  --power -.25 \
#  --bfile ${bed} \

#bash "${ldak}" \
#  --join-kins ${out_dir}

