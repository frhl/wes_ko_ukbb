# wes_ko_ukbb

Counts can be found here https://docs.google.com/spreadsheets/d/1p3wvZx7BNcMDg2BBT1wNSpVjC9gEwOloDEqeKMrcRAM/edit#gid=0

updated 17-sep-2021

## todo
* (DONE) use UKBB kinship markers to generate a sparse matrix of all chromosomes
* (DONE) use these markers to generate a sparse matrix (also filter by 95% missingness)
* (DONE) Should result in around 93k markers in total.
* Setup SAIGE pipeline (fit null)
* setup SAIGE pipeline (SPA test)
* get phenotypes (Teresa)
* extract allels that are compound hetz. Validate with known compound hetz analysis.
* run saige with with sparse matrix and fake VCF files that represent knockouts. (Notice, that we treat these files as single markers and not genesets!). Note: when running SAIGE, you can use the full file with non white british. These will be removed since GRM is created on only WB.

## Current pipeline:
* Core functions are found in utils/hail_export.py
* Functions are called in shell scripts within the /scripts folder
* remember to use either 'jupyter-hail' condaenv when working with this data



