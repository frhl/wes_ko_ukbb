
devtools::load_all("utils/modules/R/gwastools")
source("scripts/post_hoc/utils.R")
library(data.table)
library(argparse)

main <- function(args){

    chrom = args$chrom
    annotation = args$annotation
    out_prefix = args$out_prefix
    path_phenotypes = args$path_phenotypes
    path_dict = args$path_pheno_dict

    # iterate over phenotypes
    df_phenotype <- fread(path_phenotypes)
    phenotypes <- readLines(path_header)
    
    # read knockouts
    d <- read_ukb_wes_kos(annotation = annotation, chromosomes = chrom)
    d$is_cis <- d$knockout %in% "Compound heterozygote (cis)"
    d$is_chet <- d$knockout %in% "Compound heterozygote"
    d$is_hom <- d$knockout %in% "Homozygote"
    d <- d[(d$is_cis) | (d$is_chet) | (d$is_hom), ]

    out <- do.call(rbind, lapply(phenotypes, function(pheno){   
        case_eid <- df_phenotype$eid[as.logical(df_phenotype[[pheno]])]
        d_cases <- d[d$s %in% case_eid,]
        cases_unpacked <- NULL
        if (nrow(d_cases) > 0) {
            cases_packed <- gwastools::cooccur_pack_lib(d_cases)
            cases_unpacked <- gwastools::cooccur_unpack_lib(cases_packed)
            cases_unpacked$phenotype <- pheno
        }
        return(cases_unpacked)
   }))
   outfile <- paste0(out_prefix,".txt.gz")
   write(paste("writing", outfile), stderr()) 
   fwrite(out, outfile, sep = "\t")
     

}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "")
parser$add_argument("--annotation", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_phenotypes", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_pheno_dict", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_header", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


