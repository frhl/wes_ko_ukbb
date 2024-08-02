devtools::load_all("utils/modules/R/gwastools")
library(argparse)
library(data.table)

main <- function(args){


    files <- list.files(args$target_dir, full.names = TRUE, pattern = ".txt.gz")
    phenotypes <- gsub_phenotype_from_path(files)
    genes <- fread(args$path_sig_genes)

    dt <- rbindlist(lapply(files, function(f){
        pheno <- gsub_phenotype_from_path(f)
        d <- fread(f)
        subset_genes <- genes$gene[(genes$trait %in% pheno)]
        d <- d[d$MarkerID %in% subset_genes,]
        d$phenotype <- pheno
        d$prs <- grepl(pattern="locoprs", f)
        d$filepath <- f
        return(d)
    }))

    # only include hits where we condition on something
    # as most of them are other variants included in the file
    dt <- dt[!is.na(dt$p.value_c),]
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(dt, outfile, sep = "\t", quote=FALSE)

}

parser <- ArgumentParser()
parser$add_argument("--target_dir", default=NULL, help = "")
parser$add_argument("--path_sig_genes", default=NULL, help = "")
parser$add_argument("--out_prefix", default=NULL, help = "")
args <- parser$parse_args()

main(args)


