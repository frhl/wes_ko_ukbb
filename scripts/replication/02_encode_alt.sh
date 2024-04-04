# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -eu

readonly group="original"
readonly in_dir="/mnt/project/wes_ko_ukbb/data/phased/export_alt"
readonly vep_dir="/mnt/project/wes_ko_ukbb/data/vep/${group}"

for pop in "eur"; do
  out_dir="/wes_ko_ukbb/data/phased/encode_alt/${pop}/${group}/by_chrom"
  dx mkdir -p ${out_dir}
  for anno in "pLoF_damaging_missense"; do
    for pp in "0.90"; do
      for af in "05"; do
        for CHR in {1..22} X; do
          out="UKB.wes.chr${CHR}.phased.qc.${pop}.af${af}.pp${pp}.${group}.${anno}.txt.gz"
          geno="${in_dir}/UKB.wes.chr${CHR}.phased.qc.${pop}.af${af}.txt.gz"
          vep="${vep_dir}/UKB.chr${CHR}.exome_array.variants_only.vep.csqs.worst_csq_by_gene_canonical.${group}.txt.gz"
          echo "GENO = ${geno}"
          echo "VEP =  ${vep}"
          echo "OUT = ${out_dir}/${out}"
          if [[ "$(dx_is_empty $( wo_mnt_project ${geno} ))" -eq 0 ]]; then
            if [[ $(dx_file_exists "${out_dir}/${out}") -eq 0 ]]; then
              if [[ "${anno}" == "pLoF_damaging_missense" ]]; then
                dx run app-swiss-army-knife \
                  -iimage_file="/docker/call_chets.tar.gz"\
                  -icmd="
                      zcat ${vep} | grep -Ew '(pLoF)|(damaging_missense)' | cut -f3,5 | gzip > map.txt.gz &&
                      zcat ${vep} | grep -Ew '(pLoF)|(damaging_missense)' | cut -f3,10 | gzip > info.txt.gz &&
                      zcat ${geno} | awk '\$4>=${pp} || \$4==\"\" || \$4==\".\"' | gzip > geno.txt.gz &&
                      call_chets --geno geno.txt.gz --gene-map map.txt.gz --show-variants | gzip > ${out}
                      rm map.txt.gz &&
                      rm geno.txt.gz &&
                      rm info.txt.gz &&
                      echo 'xpzx ${anno}'
                  " \
                  --instance-type mem1_ssd1_v2_x8 \
                  --folder=".${out_dir}" \
                  --priority normal \
                  --name annotate_${anno}_chets_chr${CHR} -y
              else
               dx run app-swiss-army-knife \
                  -iimage_file="/docker/call_chets.tar.gz"\
                  -icmd="
                      zcat ${vep} | grep -Ew ${anno} | cut -f3,5 | gzip > map.txt.gz &&
                      zcat ${vep} | grep -Ew ${anno} | cut -f3,10 | gzip > info.txt.gz &&
                      zcat ${geno} | awk '\$4>=${pp} || \$4==\"\" || \$4==\".\"' | gzip > geno.txt.gz &&
                      call_chets --geno geno.txt.gz --gene-map map.txt.gz --show-variants | gzip > ${out}
                      rm map.txt.gz &&
                      rm geno.txt.gz &&
                      rm info.txt.gz &&
                      echo 'okzza ${anno}'
                  " \
                  --instance-type mem1_ssd1_v2_x8 \
                  --folder=".${out_dir}" \
                  --priority normal \
                  --name c${CHR}_${pop}_encode_${anno} -y
                fi
              else
                >&2 echo "${out} already exists! Skipping.."
              fi
            else
              >&2 echo "${geno} does not exists or is empty!"
            fi
          done
        done
     done
  done
done






