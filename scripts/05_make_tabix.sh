#!/usr/bin/env bash
#
#
#$ -N tabix
#$ -wd /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb
#$ -o logs/tabix.log
#$ -e logs/tabix.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 1
#$ -q short.qc@@short.hge

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/vcf_utils.sh

module load BCFtools/1.12-GCC-10.3.0

readonly in_dir="derived/knockouts/211102"

for i in "${in_dir}/"*.vcf.bgz; do
    [ -f "$i" ] || break
    echo "running make_tabix ${i}"
    make_tabix "${i}" "csi"
done

# index with csi
#make_tabix "${out_prefix}_ko.vcf.bgz" "csi"



