#!/bin/bash
set -eu

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

# we automatically upload local scripts to remote when they change
readonly remote_dir="wes_ko_ukbb/scripts"
readonly rscript_local="01_gen_carrier_matrix.R"
readonly rscript_remote="${remote_dir}/01_gene_carrier_matrix.R"
dx_update_remote ${rscript_remote} ${rscript_local}


# parameters for replication
# readonly group="replication_full"
readonly group="original"
readonly pop="eur"
readonly anno="pLoF_damaging_missense"
readonly pp="0.90" 
readonly af="05"

readonly genes_per_block="100"
readonly sample_dir="/mnt/project/wes_ko_ukbb/data/phenotypes"
readonly in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt/${pop}/${group}"
readonly out_dir="/wes_ko_ukbb/data/survival/replication/input/2402/${group}/blocks"
#readonly path_genes_to_test="/mnt/project/wes_ko_ukbb/data/replication/cox_genes_to_test.txt"
readonly path_genes_to_test="/mnt/project/wes_ko_ukbb/data/replication/2401_cox_replication_genes.txt"

dx mkdir -p ${out_dir}

for infile in "remove.wes200k" "keep.wes200k"; do
  in_file="${in_dir}/UKB.wes.merged.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.txt.gz"
  sample_file="${sample_dir}/qced_bin_matrix_${pop}.${infile}.samples"
  out_prefix="UKB.carrier_matrix.qced.${pop}.unrel.af${af}.pp${pp}.${anno}.${infile}"
  if [[ $(dx_file_exists "${out_dir}/${out_prefix}.txt") -eq 0 ]]; then
    dx run app-swiss-army-knife \
      -iimage_file="/docker/rsuite.tar.gz"\
      -icmd="
          Rscript /mnt/project/${rscript_remote} \
           --input_file ${in_file} \
           --set_genes_to_test ${path_genes_to_test} \
           --sample_file ${sample_file} \
           --out_prefix ${out_prefix} \
           --genes_per_block ${genes_per_block} \
           --n_min_ko 5
        "\
      --instance-type mem1_ssd1_v2_x8 \
      --folder=".${out_dir}" \
      --priority high \
      --name gen_carrier_matrix_${anno}_af${af}_pp${pp} -y 
  else
    >&2 echo "${out_prefix}.txt already exists."
  fi
done



