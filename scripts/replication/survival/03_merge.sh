#!/bin/bash

anc="eur"
pp="0.90"
af="05"
anno="pLoF_damaging_missense"
group="replication_full"

set -eu
in_dir="/mnt/project/wes_ko_ukbb/data/survival/replication/output/${group}"
out_dir="/wes_ko_ukbb/data/survival/replication/output"

dx mkdir -p ${out_dir}

for infile in "keep.wes200k" "remove.wes200k"; do
  in_prefix="${in_dir}/UKB.tte.chr*.unrel.af${af}.pp${pp}.${anno}.*.${infile}.txt.gz"
  out_file="UKB.tte.combined.anc${anc}.unrel.af${af}.pp${pp}.${anno}.${group}.${infile}.txt"
  dx run app-swiss-army-knife \
    -icmd="
      zcat ${in_prefix} | head -n 1 > ${out_file} &&
      for f in ${in_prefix}; do
        zcat \$f | wc -l;
        zcat \$f | tail -n+2 >> ${out_file};
      done &&
      gzip ${out_file}
      "\
    --instance-type mem1_ssd1_v2_x8 \
    --folder=".${out_dir}" \
    --priority normal \
    --name merge_coxph_${anno}_af${af}_pp${pp}_${anc} -y 
done


