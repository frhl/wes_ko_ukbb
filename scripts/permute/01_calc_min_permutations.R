#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# helper, read a hail list
clean_hail_list <- function(x, comma_sub = ', ') {
    x <- gsub('(\\[)|(\\])|(\\")','',x)
    x <- gsub('\\,',comma_sub,x)
    return(x)
}


# Read with and annotate basename
fread_with_basename <- function(f){
    d <- fread(f)
    file_sans_gz <- tools::file_path_sans_ext(basename(f))
    file_sans_txt <- tools::file_path_sans_ext(file_sans_gz)
    d$basename <- file_sans_txt
    return(d)
}

# generate a table of genes in file
tabulate_MarkerID_basename <- function(M){
    stopifnot(all(c("MarkerID",'basename') %in% colnames(M)))
    d <- as.data.frame(table(M$MarkerID, M$basename))
    colnames(d) <- c("MarkerID", "basename", "Freq")
    d <- d[d$Freq > 0,]
    d$Freq <- NULL
    return(d)
}

main <- function(args){

    print(args)
    stopifnot(dir.exists(args$spa_cts_dir))     
    stopifnot(dir.exists(args$spa_bin_dir))

    spa_cts_files <- list.files(args$spa_cts_dir, pattern = ".txt.gz", full.names = TRUE)
    spa_bin_files <- list.files(args$spa_bin_dir, pattern = ".txt.gz", full.names = TRUE)
    
    # read in all the SPA step 2 results
    spa_cts_full <- do.call(rbind, lapply(spa_cts_files, fread_with_basename))
    spa_bin_full <- do.call(rbind, lapply(spa_bin_files, fread_with_basename))

    # write out each gene involved in each comparison (avg of ~1700 genes per SPA)
    spa_full <- rbind(spa_cts_full[,c("MarkerID","basename", "p.value")], spa_bin_full[,c("MarkerID","basename", "p.value")])
    print(head(spa_full))
    out_prefix_genes <- paste0(args$out_prefix, "_genes.tsv.gz")
    write(paste("[verbose] writing", out_prefix_genes),stderr())
    fwrite(tabulate_MarkerID_basename(spa_full), out_prefix_genes, sep = '\t')
    
    # write out phenotype-genes P-value pairs
    out_prefix_p <- paste0(args$out_prefix, "_pvalue_genes.tsv.gz")
    write(paste("[verbose] writing", out_prefix_p),stderr())
    fwrite(spa_full, out_prefix_genes, sep = '\t')

    # now combine the two files and get min P-value detected for each gene
    spa_cts <- spa_cts_full[,c('MarkerID','p.value', 'CHR')]
    spa_bin <- spa_bin_full[,c('MarkerID','p.value', 'CHR')]
    spa <- rbind(spa_cts, spa_bin)
    spa$p.value <- as.numeric(spa$p.value)
    NAs <- suppressWarnings(sum(is.na(spa$p.value)))
    if (NAs > 0) write(paste("[verbose]",NAs, "NAs in SPA P.values"), stderr())
    
      # aggregate p-value by marker id and take min
    spa_chr <- spa[,c('MarkerID', 'CHR')]
    spa_chr <- spa_chr[!duplicated(spa_chr),]
    spa_aggr <- aggregate(p.value ~ MarkerID, data = spa, FUN = function(x) min(x, na.rm = TRUE))
    
    # deal with extremely low P-values (that cause Inf permutations)
    p_extreme <- spa_aggr$p.value < as.numeric(args$p_cutoff)
    spa_aggr$p.value[p_extreme] <- as.numeric(args$p_cutoff)
 
    # calculate needed number of permutations reqeuired
    spa_aggr$permut <- ceiling(1/spa_aggr$p.value)

    tsv <- do.call(rbind, lapply(1:22, function(chr){
        path <- gsub('CHR',chr, args$tsv_path)
        fread(path)
    }))

    tsv <- tsv[,c('gene_id','varid')]
    tsv <- tsv[!duplicated(tsv),]
    tsv$varid <- clean_hail_list(tsv$varid)

    aggr <- aggregate(varid ~ gene_id, data = tsv, FUN = function(x) paste(x, collapse = ', '))
    aggr$varid <- unlist(lapply(strsplit(aggr$varid, split = ', '), function(x) paste0(unique(x), collapse = ', ')))
    colnames(aggr) <- c('MarkerID','varid')

    spa_aggr <- merge(spa_aggr, aggr, all.x = TRUE) # combine SPA / tsv
    spa_aggr <- merge(spa_aggr, spa_chr, all.x = TRUE) # add chromsome
    spa_aggr <- spa_aggr[rev(order(spa_aggr$permut)),]
    
    # mark genes that are only homozygous (i.e. only one variant in the gene)
    n_variants <- unlist(lapply(strsplit(spa_aggr$varid, split = ','), length))
    spa_aggr$n_variants <- n_variants
    spa_aggr$knockout <- ifelse(n_variants > 1, 'CH', 'HOM only')

    outfile <- paste0(args$out_prefix, ".tsv.gz")
    write(paste("[verbose] writing",outfile),stderr())
    fwrite(spa_aggr, outfile, sep = '\t', scipen = 50)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--spa_cts_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--spa_bin_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--tsv_path", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--p_cutoff", default=10e-12, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

