#!/bin/bash
set -eu

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly anc="eur"
readonly anno="pLoF_damaging_missense"
readonly pp="0.90" 
readonly af="05"


# parameters for replication
#for group in "replication_full" "replication"; do
for group in "original"; do
  in_dir="/wes_ko_ukbb/data/survival/replication/input/2402/${group}/blocks"
  out_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/survival/from_rap/240226"
  mkdir -p ${out_dir}
  for infile in "remove.wes200k" "keep.wes200k"; do
    the_file="UKB.carrier_matrix.qced.${anc}.unrel.af${af}.pp${pp}.${anno}.${infile}.b1of1.txt.gz"
    out_file="UKB.carrier_matrix.qced.${anc}.unrel.af${af}.pp${pp}.${anno}.${infile}.${group}.b1of1.txt.gz"
    file_to_fetch="${in_dir}/${the_file}"
    dx download ${file_to_fetch}
    mv ${the_file} ${out_dir}/${out_file} 
  done
done


