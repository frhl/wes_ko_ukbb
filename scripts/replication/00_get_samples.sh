# get replication samples based on 200k exome release

readonly rap_qced_samples="/Barney/qc/09_0_final_sample_qc/09_final_qc.keep.BRaVa.sample_list"
readonly bmrc_200k_samples="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/phenotypes/samples_176k.txt"

readonly pop="eur"
readonly pop_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/samples"
readonly pop_file="${pop_dir}/UKB.chr21.samples.${pop}.txt"

readonly out_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/samples"
readonly out_prefix="${out_dir}/UKB.qced_samples.wes450k.wes200k"
readonly qced_samples="${out_dir}/UKB.qced_samples.wes450k.txt"

set -x

# Donwload QCed samples
dx cat ${rap_qced_samples} | sort > ${qced_samples}

# Extract samples exclusive to rap_qced_samples
comm -23 <(cat ${qced_samples} | sort ) <(cat ${bmrc_200k_samples} | sort) > ${out_prefix}.remove.txt
comm -12 <(cat ${out_prefix}.remove.txt | sort) <(cat ${pop_file} | sort) > ${out_prefix}.${pop}.remove.txt

# Extract overlapping samples
comm -12 <(dx cat ${rap_qced_samples} | sort ) <(cat ${bmrc_200k_samples} | sort) > ${out_prefix}.intersect.txt


