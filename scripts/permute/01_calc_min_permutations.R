#!/usr/bin/env Rscript

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

    print(args)
    #stopifnot(dir.exists(args$spa_cts_dir))
    stopifnot(dir.exists(args$spa_bin_dir))

    #spa_cts_files <- list.files(args$spa_cts_dir, pattern = ".txt.gz", full.names = TRUE)
    spa_bin_files <- list.files(args$spa_bin_dir, pattern = ".txt.gz", full.names = TRUE)

    # read in all the SPA step 2 results
    #spa_cts_full <- rbindlist(lapply(spa_cts_files, fread_with_basename), fill = TRUE)
    spa_bin_full <- rbindlist(lapply(spa_bin_files, fread_with_basename), fill = TRUE)
    write("Loaded all cts/binary files..", stdout())

    #print(head(spa_cts_full, n=2))
    print(head(spa_bin_full, n=2))
    
        # combine and subset
    keep <- c("CHR","MarkerID","basename", "p.value", "Tstat", "p.value_c", "Tstat_c")
    #spa_full <- rbindlist(list(spa_cts_full, spa_bin_full))
    spa_full <- spa_bin_full
    spa_full <- spa_full[, colnames(spa_full) %in% keep, with = FALSE]
    spa_full <- spa_full[,..keep]
    
    # what is the format of the input file?
    basename_prefix = args$basename_prefix
    basename_annotation = "(pLoF_damaging_missense)|(pLoF)|(damaging_missense)"
    basename_locoprs = "locoprs"

    # perform string operations to get pheotype
    stopifnot(all(grepl(basename_prefix, spa_full$basename)))
    stopifnot(all(grepl(basename_annotation, spa_full$basename)))
    basename_minimal <- gsub(basename_prefix, "", spa_full$basename)
    basename_minimal <- gsub(basename_annotation, "", basename_minimal)
    basename_minimal <- gsub(basename_locoprs, "", basename_minimal)
    spa_full$phenotype <- gsub("(^\\_*)|(\\_*$)", "", basename_minimal)
    spa_full$annotation <- stringr::str_extract_all(spa_full$basename, basename_annotation)
    spa_full$prs <- ifelse(grepl("locoprs", spa_full$basename), "locoprs", NA)
    spa_full$basename <- NULL

    n <- nrow(spa_full) 

    if (args$use_cond_p) {
        
       # use p-value from conditional values if applicable 
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

    # discard markers from files that are significant
    spa_full <- spa_full[grepl("ENSG", spa_full$MarkerID),]

    # lookup hashes and AC
    genes <- unique(spa_full$MarkerID) # ~ 1800 genes
    phenos <- unique(spa_full$phenotype) # ~ 45 phenotypes

    # get AC for relevant phenotypes and 
    AUTOSOMES <- 1:22 
    lst_ac <- lapply(AUTOSOMES, function(chr){
        path <- gsub("CHR",chr,args$ac_path)
        if (!file.exists(path)) stop(paste(path, "does not exist"))
        d <- fread(path)
        d$chr <- NULL
        d$ref <- NULL
        d$alt <- NULL
        d$pos <- NULL
        return(d)
    })
    
    # Get hash 
    lst_hash <- lapply(AUTOSOMES, function(chr){
        path <- gsub("CHR",chr,args$hash_path)
        if (!file.exists(path)) stop(paste(path, "does not exist"))
        d <- fread(path)
        d$chr <- NULL
        d$ref <- NULL
        d$alt <- NULL
        d$pos <- NULL
        return(d)
    })

    # merge allele count onto main file
    d_ac <- do.call(rbind, lst_ac)
    d_hash <- do.call(rbind, lst_hash)
    melted_ac <- data.table::melt(d_ac, id.var = "id")
    melted_hash <- data.table::melt(d_hash, id.var = "id")
    colnames(melted_ac) <- c("MarkerID", "phenotype", "AC")
    colnames(melted_hash) <- c("MarkerID", "phenotype", "HASH")
    spa_full <- merge(spa_full, melted_hash, all.x = TRUE)
    spa_full <- merge(spa_full, melted_ac, all.x = TRUE)

    # format to avoid scientific notation 
    out_prefix_true_p_detailed <- paste0(args$out_prefix, "_true_p_detailed.tsv.gz")
    fwrite(spa_full, out_prefix_true_p_detailed, sep = '\t', quote = FALSE, na = NA)

    # remove old columns
    #spa_full$p.value_c <- NULL
    #spa_full$Tstat_c <- NULL
    spa_full$cond <- NULL

        # the filtered version
    #spa_full$pvalue <- format(spa_full$pvalue, scientific = FALSE)
    #spa_full$tstat <- format(spa_full$tstat, scientific = FALSE)
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
parser$add_argument("--spa_cts_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--spa_bin_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--ac_path", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--hash_path", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--basename_prefix", default="ukb_eur_wes_200k_maf0to5e-2", required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--tsv_path", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--p_cutoff", default=NULL, required = FALSE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--use_cond_p", default=FALSE, action = "store_true", help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

