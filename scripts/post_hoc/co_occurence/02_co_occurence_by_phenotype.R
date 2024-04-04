
devtools::load_all("utils/modules/R/gwastools")
source("scripts/post_hoc/utils.R")
library(data.table)
library(argparse)

main <- function(args){

    chrom = args$chrom
    annotation = args$annotation
    out_prefix = args$out_prefix
    path_phenotypes = args$path_phenotypes
    path_header = args$path_header
    convert_tte_to_bool = args$convert_tte_to_bool

    # iterate over phenotypes
    df_phenotype <- fread(path_phenotypes)
    phenotypes <- readLines(path_header)

    # for testing, just take first 30
    write("Warning! We are only doing first 10 phenos for testing. Edit code to change!", stderr())
    phenotypes <- phenotypes[1:10]

    # convert numerics to bools for time-to-event phenotypes
    if (convert_tte_to_bool){
        eids <- df_phenotype$eid
        cols <- which(colnames(df_phenotype) %in% phenotypes)
        df_phenotype[cols] <- lapply(df_phenotype[cols], function(x) !is.na(suppressWarnings(as.numeric(x))))        
        # TTE controls are NAs
        df_phenotype[!is.na(df_phenotype)] <- TRUE
        df_phenotype[is.na(df_phenotype)] <- FALSE 
        df_phenotype$eid <- eids
    }

    # read knockouts
    d <- read_ukb_wes_kos(annotation = annotation, chromosomes = chrom)
    d$is_cis <- d$knockout %in% "Compound heterozygote (cis)"
    d$is_chet <- d$knockout %in% "Compound heterozygote"
    d$is_hom <- d$knockout %in% "Homozygote"
    d$is_het <- d$knockout %in% "Heterozygote"
    d <- d[(d$is_cis) | (d$is_chet) | (d$is_hom) | (d$is_het), ]

    # read over cases that are knockouts
    out <- do.call(rbind, lapply(phenotypes, function(pheno){   
        # deal with cases
        print(pheno)
        case_eid <- df_phenotype$eid[as.logical(df_phenotype[[pheno]])]
        print(paste("cases:", length(case_eid)))
        d_cases <- d[d$s %in% case_eid,]
        cases_unpacked <- NULL
        if (nrow(d_cases) > 0) {
            cases_packed <- gwastools::cooccur_pack_lib(d_cases)
            cases_unpacked <- gwastools::cooccur_unpack_lib(cases_packed)
            cases_unpacked$phenotype <- pheno
            cases_unpacked$type <- "cases"
        }
        # deal with controls and cases (i.e. full knockout count)
        control_eid <- df_phenotype$eid[!is.na(df_phenotype[[pheno]])]
        print(paste("controls:", length(control_eid)))
        case_control_eid <- c(case_eid, control_eid)
        d_cases_controls <- d[d$s %in% case_control_eid,]
        cases_controls_unpacked <- NULL
        if (nrow(d_cases_controls) > 0){
            cases_controls_packed <- gwastools::cooccur_pack_lib(d_cases_controls)
            cases_controls_unpacked <- gwastools::cooccur_unpack_lib(cases_controls_packed)
            cases_controls_unpacked$phenotype <- pheno
            cases_controls_unpacked$type <- "all"
        }

        combined <- rbind(cases_unpacked, cases_controls_unpacked)
        return(combined)
   }))
   outfile <- paste0(out_prefix,".txt.gz")
   write(paste("writing", outfile), stderr()) 
   fwrite(out, outfile, sep = "\t")
}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "")
parser$add_argument("--annotation", default=NULL, required = TRUE, help = "")
parser$add_argument("--convert_tte_to_bool", default=FALSE, action="store_true", help = "")
parser$add_argument("--path_phenotypes", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_header", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


