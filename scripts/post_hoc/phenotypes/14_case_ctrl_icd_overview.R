
library(data.table)
library(argparse)

main <- function(args){
        
    stopifnot(file.exists(args$ko_by_phenotype1))
    stopifnot(file.exists(args$ko_by_phenotype2))
    stopifnot(file.exists(args$ko_by_phenotype3))
    stopifnot(file.exists(args$icd_path))


    # read pLoF_damaging_missense
    d1 <- fread(args$ko_by_phenotype1)
    colnames(d1)[1:5] <- c("unix_code","cis","chet","het","hom")
    d1 <- d1[d1$controls_nfe == 0,]
    d1$controls_nfe <- NULL 
    d1$cases_nfe <- NULL 
    colnames(d1)[-1] <- paste0("pLoF_damaging_missense_",colnames(d1)[-1])

    # read pLoFs
    d2 <- fread(args$ko_by_phenotype2)
    colnames(d2)[1:5] <- c("unix_code","cis","chet","het","hom")
    d2 <- d2[d2$controls_nfe == 0,]
    d2$controls_nfe <- NULL 
    d2$cases_nfe <- NULL 
    colnames(d2)[-1] <- paste0("pLoF_",colnames(d2)[-1])

    # read damaging_missense
    d3 <- fread(args$ko_by_phenotype3)
    colnames(d3)[1:5] <- c("unix_code","cis","chet","het","hom")
    d3 <- d3[d3$controls_nfe == 0,]
    d3$controls_nfe <- NULL
    d3$cases_nfe <- NULL 
    colnames(d3)[-1] <- paste0("damaging_missense_",colnames(d3)[-1])

    # combine all of them
    mrg <- merge(d2, d3, by = c("unix_code"))
    mrg <- merge(mrg , d1, by = c("unix_code"))
    mrg <- mrg[!duplicated(mrg),]

    # combine them and rename
    stopifnot(all(mrg$pLoF_damaging_missense_cases_all == mrg$pLoF_phenos_cases_all))
    stopifnot(all(mrg$pLoF_damaging_missense_cases_all == mrg$damaging_missense_phenos_cases_all))
    stopifnot(all(mrg$pLoF_damaging_missense_controls_all == mrg$pLoF_phenos_controls_all))
    stopifnot(all(mrg$pLoF_damaging_missense_controls_all == mrg$damaging_missense_phenos_controls_all))
    mrg$damaging_missense_phenos_na <- NULL
    mrg$damaging_missense_cases_all <- NULL
    mrg$damaging_missense_controls_all <- NULL
    mrg$pLoF_phenos_na <- NULL
    mrg$pLoF_cases_all <- NULL
    mrg$pLoF_controls_all <- NULL

    # rename since they are all the same
    colnames(mrg)[colnames(mrg)=="pLoF_damaging_missense_phenos_na"] <- "undefined"
    colnames(mrg)[colnames(mrg)=="pLoF_damaging_missense_cases_all"] <- "cases"
    colnames(mrg)[colnames(mrg)=="pLoF_damaging_missense_controls_all"] <- "ctrls"

    # merge with ICD codes
    icd <- fread(args$icd_path)
    mrg <- merge(icd, mrg, by = "unix_code", all.x = TRUE)
    mrg <- mrg[rev(order(mrg$ICD_chapter)),]
        
    
    # write file containing all annotations 
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(mrg, outfile, sep = "\t", na="NA")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--ko_by_phenotype1", default=NULL, required = TRUE, help = "path to knockout by phenotypes (e.g. pLoFs)")
parser$add_argument("--ko_by_phenotype2", default=NULL, required = TRUE, help = "--")
parser$add_argument("--ko_by_phenotype3", default=NULL, required = TRUE, help = "--")
parser$add_argument("--icd_path", default=NULL, required = TRUE, help = "--")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


