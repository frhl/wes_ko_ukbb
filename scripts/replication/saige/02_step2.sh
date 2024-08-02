dx build -f saige-universal-step-2

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset

readonly pheno_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/phenotypes"
readonly step0_dir="/wes_ko_ukbb/data/saige/step0"
readonly step1_dir="/wes_ko_ukbb/data/saige/replication/step1"
readonly test_type="variant"

readonly pop="eur"
readonly group="original"
readonly out_dir="/wes_ko_ukbb/data/saige/replication/polished/step2_qced/${group}"
readonly exome_dir="/wes_ko_ukbb/data/phased/encode_alt/${pop}/${group}"

readonly af="05"
readonly pp="0.90"
readonly mode="recessive"
readonly anno="pLoF_damaging_missense"
   
dx mkdir -p ${out_dir}
readonly bin_phenos="${pheno_dir}/bin_matrix_eur_header_replication_phenos.txt"
for pheno in $(cat $bin_phenos); do
  for infile in "keep.wes200k" "remove.wes200k"; do 
    #step1_prefix="${pheno}.${infile}"
    step1_prefix="${pheno}.qced.${infile}"
    exome_prefix="${exome_dir}/UKB.wes.merged.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.${mode}.auto"
    out_prefix="UKB.wes.${pheno}.merged.phased.qc.${infile}.${pop}.af${af}.pp${pp}.${group}.${anno}.${mode}.auto"
    if [[ $(dx_file_exists "${out_dir}/${out_prefix}.txt.gz") -eq "0" ]]; then
       dx run saige-universal-step-2 \
                -i chrom="1" \
                -i output_prefix="${out_prefix}" \
                -i model_file="${step1_dir}/${step1_prefix}.rda" \
                -i variance_ratio="${step1_dir}/${step1_prefix}.varianceRatio.txt" \
                -i test_type=$test_type \
                -i exome_bed="${exome_prefix}.bed" \
                -i exome_bim="${exome_prefix}.bim" \
                -i exome_fam="${exome_prefix}.fam" \
                -i GRM="${step0_dir}/ukb_array_400k_eur_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx" \
                -i GRM_samples="${step0_dir}/ukb_array_400k_eur_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
                --instance-type "mem3_ssd1_v2_x8" --priority normal --destination "${out_dir}" -y --name "s2.${pheno}.${infile}"
     else 
        >&2 echo "${out_prefix}.txt.gz already exists. Skipping.."
     fi
  done
done

