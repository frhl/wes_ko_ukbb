# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly group="original"
for pop in "eur"; do
  in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt/${pop}/${group}"
  out_dir="/wes_ko_ukbb/data/phased/encode_alt/${pop}/${group}"
  for mode in "additive" "recessive"; do
    for anno in "pLoF_damaging_missense"; do
      for pp in "0.90"; do
        for af in "05"; do
          in_vcf_auto="${in_dir}/UKB.wes.merged.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.${mode}.auto.vcf.gz"
          in_vcf_chrx="${out_dir}/UKB.wes.merged.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.${mode}.chrx.vcf.gz"
          out_prefix_auto="UKB.wes.merged.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.${mode}.auto"
          out_prefix_chrx="UKB.wes.merged.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.${mode}.chrx"
          if [[ $(dx_file_exists "${out_dir}/${out_prefix_auto}.bed") -eq "0" ]]; then 
            test_vcf="$(wo_mnt_project ${in_vcf_auto})"
            if [[ $(dx_file_exists ${test_vcf}) -eq "1" ]]; then
              dx run app-swiss-army-knife \
                -icmd="
                  plink2 --vcf ${in_vcf_auto} dosage=DS --make-bed --out ${out_prefix_auto} && echo 'ok' 
                  "\
                --instance-type mem1_ssd1_v2_x8 \
                --folder=".${out_dir}" \
                --priority normal \
                --name vcf_to_plink_${mode}_${anno} -y
            else
               >&2 echo "Error! '${in_vcf_auto}' does not exists."
            fi
          fi
        done
      done
    done
  done
done

