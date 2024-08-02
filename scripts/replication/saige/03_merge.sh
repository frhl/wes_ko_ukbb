# cohort dataset
source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

# keep track of updates to rscript
remote_dir="wes_ko_ukbb/scripts"
rscript_local="03_merge.R"
rscript_remote="${remote_dir}/03_merge.R"
dx_update_remote ${rscript_remote} ${rscript_local}

set -o errexit
set -o nounset

gene_map="/mnt/project/genesets/220524_hgnc_ensg_enst_chr_pos.txt.gz"

modality="binary"
#group="replication_full"
group="original"
out_dir="/wes_ko_ukbb/data/saige/replication/merged"
#in_dir="/mnt/project/wes_ko_ukbb/data/saige/replication/step2_wes200k/${group}"
in_dir="/mnt/project/wes_ko_ukbb/data/saige/replication/polished/step2_qced/${group}"
#input_pattern="UKB.*.txt.gz"
input_pattern="UKB.*.txt.gz"

out_prefix="UKB.saige.${modality}.${group}.240226.replication.original.final"


dx mkdir -p ${out_dir}
dx run app-swiss-army-knife \
  -iimage_file="/docker/rsuite.tar.gz"\
  -icmd="Rscript /mnt/project/${rscript_remote} \
     --out_prefix '${out_prefix}' \
     --input_dir '${in_dir}' \
     --pattern '${input_pattern}' \
     --gene_map_path '${gene_map}' &&
     echo 'ok'
  " \
  --instance-type mem2_ssd1_v2_x32 \
  --folder=".${out_dir}" \
  --priority normal \
  --name merge_plink -y               



