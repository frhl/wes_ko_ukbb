#!/bin/bash

main() {
    
    dx-mount-all-inputs --except GRM --except genotype_bed --except genotype_bim --except genotype_fam 

    dx download "$genotype_bed" -o genotype.bed
    dx download "$genotype_bim" -o genotype.bim
    dx download "$genotype_fam" -o genotype.fam
    dx download "$GRM" -o GRM

    git clone https://github.com/BRaVa-genetics/universal-saige

    mv universal-saige/* .

    bash download_resources.sh --saige-image

    mkdir out

    #--invNormalize TRUE \
    ls -al

    bash 01_step1_fitNULLGLMM.sh \
        -t $trait_type \
        --genotypePlink "genotype" \
        --phenoFile in/pheno_list/* \
        --phenoCol $pheno \
        --covarColList $covariates \
        --categCovarColList "${categorical_covariates}" \
        --sampleIDs in/sample_ids/* \
        --sampleIDCol "eid" \
        --outputPrefix ${output_prefix} \
        --isSingularity false \
        --sparseGRM GRM \
        --sparseGRMID in/GRM_samples/*

    mkdir -p ~/out/model_file/
    mkdir -p ~/out/variance_ratio/

    mv ${output_prefix}.rda ""$HOME"/out/model_file/"
    mv ${output_prefix}.varianceRatio.txt ""$HOME"/out/variance_ratio/"

    dx-upload-all-outputs
}
