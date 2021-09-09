# wes_ko_ukbb

updated 9-sep-2021

## todo
* use UKBB kinship markers to generate a sparse matrix of all chromosomes
* use these markers to generate a sparse matrix (also filter by 95% missingness)
* Should result in around 93k markers in total.
* run saige with with sparse matrix and fake VCF files that represent knockouts. (Notice, that we treat these files as single markers and not genesets!).


## Current pipeline:
* Core functions are found in utils/hail_export.py
* Functions are called in shell scripts within the /scripts folder
* remember to use either 'jupyter-hail' condaenv when working with this data



