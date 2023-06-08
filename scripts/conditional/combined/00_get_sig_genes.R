devtools::load_all("utils/modules/R/gwastools")
library(argparse)
library(data.table)


main <- function(args){

    args$p_cutoff <- as.numeric(args$p_cutoff)
    files <- gwastools::list_files_saige(args$cond_step, prs=args$prs)
    print(files)
    phenos <- gwastools::gsub_phenotype_from_path(files)
    files <- files[phenos %in% gwastools::get_phenos_tested()]

    d <- rbindlist(lapply(files, function(f){
        trait <- basename(f)
        trait <- stringr::str_extract(trait, "200k_.+pLoF_damaging_missense")
        trait <- gsub("200k_", "", trait)
        trait <- gsub("_pLoF_damaging_missense", "", trait)
        d <- fread(f)
        stopifnot("MarkerID" %in% colnames(d))
        stopifnot("N_case_hom" %in% colnames(d))
        d <- d[grepl("ENSG", d$MarkerID),]
        d <- d[(d$N_case_hom > 0),]
        if ("p.value_c" %in% colnames(d)) {
            d$p.value <- d$p.value_c 
            d$cond <- TRUE
        } else {
            d$cond <- FALSE
        }
        d <- d[d$p.value < args$p_cutoff,]
        if (nrow(d)){
            out <- data.table(
                trait = trait,
                chromosome = d$CHR,
                gene = d$MarkerID,
                pvalues = d$p.value,
                cond = d$cond
            )
            return(out)
        }
    }))
    # write files
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(d, outfile, sep = "\t", quote=FALSE)

}

parser <- ArgumentParser()
parser$add_argument("--cond_step", default=NULL, help = "")
parser$add_argument("--prs", default=NULL, help = "")
parser$add_argument("--p_cutoff", default=5e-7, help = "")
parser$add_argument("--out_prefix", default=NULL, help = "")
args <- parser$parse_args()

main(args)


