#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
library(argparse)
library(data.table)

# Read with and annotate basename
fread_with_basename <- function(f){
    d <- fread(f)
    if (! "V1" %in% colnames(d)){
        file_sans_gz <- tools::file_path_sans_ext(basename(f))
        file_sans_txt <- tools::file_path_sans_ext(file_sans_gz)
        d$basename <- file_sans_txt
        return(d)
    } else {
        warning(paste(f, "had NA's in colname! Skipping.."))
        return(NULL)
    }
}

main <- function(args){

    # load files
    files <- gwastools::list_files_saige(cond = args$cond, prs = args$prs)
    spa_full <- rbindlist(lapply(files, fread_with_basename), fill = TRUE)

    # what is the format of the input file?
    basename_prefix = args$basename_prefix
    basename_annotation = "(pLoF_damaging_missense)|(pLoF)|(damaging_missense)"
    basename_locoprs = "locoprs"

    # perform string operations to get pheotype
    diagnosis <- stringr::str_extract(spa_full$basename, "200k_([A-Z]|[a-z]|[0-9]|\\_)+_pLoF")
    spa_full$phenotype <- gsub("(200k_)|(_pLoF)","", diagnosis)
    spa_full$annotation <- "pLoF_damaging_missense"    
    spa_full$prs <- ifelse(grepl("locoprs", spa_full$basename), "locoprs", NA)

    # append to column if it does not exists
    if (!"p.value_c" %in% colnames(spa_full)) spa_full$p.value_c <- NA
    if (!"Tstat_c" %in% colnames(spa_full)) spa_full$Tstat_c <- NA

    if (args$use_cond_p) {
       bool_c <- !is.na(spa_full$p.value_c)
       spa_full$pvalue <- spa_full$p.value
       spa_full$tstat <- spa_full$Tstat
       spa_full$cond <- bool_c
       spa_full$pvalue[bool_c] <- spa_full$p.value_c[bool_c] 
       spa_full$tstat[bool_c] <- spa_full$Tstat_c[bool_c]
    }    

    # if P-value cutoff, then discard genes
    if (!is.null(args$p_cutoff)){
        p_cutoff <- as.numeric(args$p_cutoff)
        stopifnot(p_cutoff < 1)
        bool_discard <- spa_full$pvalue > p_cutoff
        n_discard <- sum(bool_discard)
        write(paste("Removing", n_discard, " genes that do not pass P-value threshold"), stdout())
        spa_full <- spa_full[!bool_discard,] 
    }
    
    keep <- c("CHR","MarkerID","basename", "p.value", "Tstat", "p.value_c", "Tstat_c")
    spa_full <- spa_full[,..keep]

    # discard markers from files that are significant
    spa_full <- spa_full[grepl("ENSG", spa_full$MarkerID),]

    # lookup hashes and AC
    genes <- unique(spa_full$MarkerID) # ~ 1800 genes
    phenos <- unique(spa_full$phenotype) # ~ 45 phenotypes

    # format to avoid scientific notation 
    out_prefix_true_p_detailed <- paste0(args$out_prefix, "_true_p_detailed.tsv.gz")
    fwrite(spa_full, out_prefix_true_p_detailed, sep = '\t', quote = FALSE, na = NA)

    # remove old columns
    spa_full$cond <- NULL

    # the filtered version
    out_prefix_true_p <- paste0(args$out_prefix, "_true_p.tsv.gz")
    fwrite(spa_full, out_prefix_true_p, sep = '\t', quote = FALSE, na = NA)
    write(paste("Note: wrote", out_prefix_true_p),stderr())
  
    # write out genes and chromosomes used 
    spa_genes <- spa_full[,c("MarkerID","CHR")]
    spa_genes <- spa_genes[!duplicated(spa_genes),]
    out_prefix_genes <- paste0(args$out_prefix, "_genes.tsv.gz")
    fwrite(spa_genes, out_prefix_genes, sep = '\t')
    write(paste("Note: wrote", out_prefix_genes),stderr())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--cond", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--prs", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--p_cutoff", default=NULL, required = FALSE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--use_cond_p", default=FALSE, action = "store_true", help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

