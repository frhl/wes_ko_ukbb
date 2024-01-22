#!/usr/bin/env Rscript

devtools::load_all('utils/modules/R/phasingtools')
library(argparse)
library(data.table)
library(Hmisc)
library(dplyr)

# a function to get binomial CIs
bin_conf <- function(counts, col){
    counts_ci <- do.call(rbind, lapply(1:nrow(counts), function(i) Hmisc::binconf(counts[[col]][i], counts$count[i])))
    colnames(counts_ci) <- tolower(colnames(counts_ci))
    counts <- cbind(counts, counts_ci)
    return(counts)
}


main <- function(args){

    sites <- args$sites
    samples <- args$samples
    input_path <- args$input_path
    out_prefix <- args$out_prefix
    summary_type <- args$summary_type   
 
    exome_vars <- fread(sites)
    qc_samples <- readLines(samples)

    dt_summary <- list()
    for (chr in seq(1,22))
    {
        fpath <- gsub("CHR", chr, input_path) 
        cat(paste('chromosome', chr, "-", fpath))
        stopifnot(file.exists(fpath))
        chrom <- paste0("chr", chr)

        dt <- fread(fpath, key=c("locus", "alleles"))

        # deal with AC looking like this [156]
        dt$AC <- as.numeric(gsub("(\\[)|(\\])", "", dt$AC))
        dt$MAC <- dt$AC #min(dt$AC, (200000*2)-dt$AC)

        
        dt <- merge(dt, exome_vars)
        dt <- dt[dt$s %in% qc_samples,]

        rb_sizes <- dt %>% group_by(s, PS_rb) %>% summarise(count=n())
        pairs <- data.table(rb_sizes %>% filter(count == 2))
        others <- data.table(rb_sizes %>% filter(count > 2))

        setkeyv(dt, c("s", "PS_rb"))
        setkeyv(pairs, c("s", "PS_rb"))
        setkeyv(others, c("s", "PS_rb"))

        dt_pairs <- merge(dt, pairs)
        dt_others <- merge(dt, others)

        # Agreement of dt_pairs
        dt_pairs$PP[is.na(dt_pairs$PP)] <- 1
        dt_summary[[chrom]] <- dt_pairs %>% group_by(s, PS_rb) %>% 
            summarise(
            match = (sum(GT == GT_rb) != 1),
            min_MAC = min(MAC),
            max_MAC = max(MAC),
            min_PP = min(PP),
            max_PP = max(PP),
            chrom = chrom
            )
    }

    dt_summary <- rbindlist(dt_summary) 
    outfile <- paste0(out_prefix,".raw.txt.gz")
    fwrite(dt_summary, outfile, sep = '\t')

    # combinme
    if (summary_type=="shapeit"){
    setDT(dt_summary)[, MAC_bin := as.integer(cut(min_MAC, breaks = c(1, 2, 6, 11, 21, 51, 101, 201, 501, 1000, 2000, 5000, 10000, (max(dt_summary$min_MAC)+1)), right=FALSE))][,
             MAC_group := c(
             'Singleton',
             '2-5',
             '6-10',
             '11-20',
             '21-50',
             '51-100',
             '101-200',
             '201-500',
             '501-1000',
             '1k-2k',
             '2k-5k',
             '5k-10k',
             '10k+')[MAC_bin]]
    } else if (summary_type=="long") {
        setDT(dt_summary)[, MAC_bin := as.integer(cut(min_MAC, breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 21, 51, 101, 201, 501, 1000, 2000, 5000, 10000, (max(dt_summary$min_MAC)+1)), right=FALSE))][,
             MAC_group := c(
             'Singleton',
             '2',
             '3',
             '4',
             '5',
             '6',
             '7',
             '8',
             '9',
             '10',
             '11-15',
             '16-20',
             '21-50',
             '51-100',
             '101-200',
             '201-500',
             '501-1000',
             '1k-2k',
             '2k-5k',
             '5k-10k',
             '10k+')[MAC_bin]]   
    }

    dt_summary_list <- list()
    for (min_pp in c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99)) {
        dt_summary_list[[as.character(min_pp)]] <- dt_summary %>% filter(min_PP > min_pp) %>% 
            group_by(MAC_bin, MAC_group) %>% summarise(
                min_PP = min_pp,
                #mean=1-mean(as.numeric(match)),
                matches=sum(as.numeric(match)),
                count=n()
            )
    }

    # combine by PP and MAC bins
    dt_summary_pp <- rbindlist(dt_summary_list)
    dt_summary_pp$errors <- dt_summary_pp$count - dt_summary$matches
    dt_summary_pp <- bin_conf(dt_summary_pp, col = "matches")
    
    # write final file
    outfile <- paste0(out_prefix,".pp.txt.gz")
    fwrite(dt_summary_pp, outfile, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, required = TRUE, help = "path to file with PS_rb, GT_rb and GT sites.")
parser$add_argument("--sites", default=NULL, required = TRUE, help = "path to quality controlled variant sites")
parser$add_argument("--samples", default=NULL, required = TRUE, help = "path to samples to subset to")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "out prefix for files.")
parser$add_argument("--summary_type", default="shapeit", required = FALSE, help = "out prefix for files.")
args <- parser$parse_args()

main(args)

