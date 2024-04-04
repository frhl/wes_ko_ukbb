dx build -f saige-universal-step-1

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly pheno_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/phenotypes"

readonly plink_dir="/wes_ko_ukbb/data/saige/step0"
readonly pheno_mat_dir="/wes_ko_ukbb/data/phenotypes"
readonly sample_dir="/wes_ko_ukbb/data/phenotypes"
readonly out_dir="/wes_ko_ukbb/data/saige/replication/step1"

set -o errexit
set -o nounset

dx mkdir -p ${out_dir}

readonly bin_phenos="${pheno_dir}/bin_matrix_eur_header_replication_phenos.txt"
for pheno in $(cat $bin_phenos ); do

  for infile in "keep.wes200k" "remove.wes200k"; do
      
    pheno_matrix_prefix="${pheno_mat_dir}/qced_bin_matrix_eur.${infile}"
    pheno_matrix="${pheno_matrix_prefix}.txt"
    sample_id_path="${pheno_matrix_prefix}.samples"
    out_prefix="${pheno}.qced.${infile}"

     if [[ $(dx_file_exists "${out_dir}/${out_prefix}.rda") -eq "0" ]]; then
       echo "Starting ${out_dir}/${out_prefix}.rda.."
       dx run saige-universal-step-1 \
                -i output_prefix="out/${out_prefix}" \
                -i sample_ids="${sample_id_path}" \
                -i genotype_bed="${plink_dir}/ukb_array_400k_eur.plink_for_var_ratio.bed" \
                -i genotype_bim="${plink_dir}/ukb_array_400k_eur.plink_for_var_ratio.bim" \
                -i genotype_fam="${plink_dir}/ukb_array_400k_eur.plink_for_var_ratio.fam" \
                -i pheno_list="${pheno_matrix}" \
                -i pheno="$pheno" \
                -i GRM="${plink_dir}/ukb_array_400k_eur_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx" \
                -i GRM_samples="${plink_dir}/ukb_array_400k_eur_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
                -i covariates="sex,age,age2,age_sex,ukbb.centre,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10" \
                -i categorical_covariates="sex" \
                -i trait_type="binary" \
                --instance-type "mem3_ssd1_v2_x4" --priority low --destination ".${out_dir}" -y --name "s1.${pheno}.${infile}"

      else
        >&2 echo "${out_prefix}.rda already exists. Skipping.."

      fi  
  done

done

