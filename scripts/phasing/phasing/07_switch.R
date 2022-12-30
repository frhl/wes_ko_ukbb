#!/usr/bin/env Rscript

devtools::load_all('utils/modules/R/phasingtools')
library(argparse)
library(RColorBrewer)
library(stringr)
library(ggsci)

fread_phased_sites <- function(file, ...){

    # get details about chunks
    bname <- basename(file)
    method <- unlist(strsplit(bname, split = '_'))[1]

    # append to data.table
    d <- fread(file, ...)
    cnames <- colnames(d)
    d$AC <- d$info.AC
    d$AN <- d$info.AN
    d$info.AC <- NULL
    d$info.AN <- NULL

    # check that columns are OK
    d$AN_m_AC <- as.numeric(d$AN - d$AC)
    d$MAC <- as.numeric(apply(d[,c("AC","AN_m_AC")], 1, min))
    d$locus <- paste0(d$CHR,":",d$POS)
    d$method <- method
    return(d)

}

# aggregate switch errors by chromosome and minor allele frequency bins
aggregate_by_chrom <- function(files, variants){
    lst <- lapply(files, function(file){
        d <- fread_phased_sites(file)
        d$wes_variant <- d$locus %in% variants$locus
        counts <- aggregate(switches ~ wes_variant + CHR, data = d, FUN = sum)
        tested <- aggregate(switches ~ wes_variant + CHR, data = d, FUN = length)
        counts <- data.table(counts, tested = tested$switches)
        return(counts)
    })
    return(lst)
}

# aggregate by MAC by CHROM
aggregate_by_mac_chrom <- function(files, mac_bins, mac_bins_labels, variants){
    stopifnot((length(mac_bins)-1) == length(mac_bins_labels))
    stopifnot(all(mac_bins>=0))
    lst <- lapply(files, function(file){
        d <- fread_phased_sites(file)
        d$wes_variant <- d$locus %in% variants$locus
        d$mac_bin <- cut(d$MAC, breaks = mac_bins, labels = mac_bins_labels)
        counts <- aggregate(switches ~ wes_variant + mac_bin + CHR, data = d, FUN = sum)
        tested <- aggregate(switches ~ wes_variant + mac_bin + CHR, data = d, FUN = length)
        counts <- data.table(counts, tested = tested$switches)
        return(counts)
    })
    return(lst)
}

# aggregate by across chromsomes
aggregate_by_mac <- function(files, mac_bins, mac_bins_labels, variants){
    stopifnot((length(mac_bins)-1) == length(mac_bins_labels))
    stopifnot(all(mac_bins>=0))
    dt <- do.call(rbind, lapply(files, function(f){
        d <- fread_phased_sites(f)
        d$wes_variant <- d$locus %in% variants$locus
        d$mac_bin <- cut(d$MAC, breaks = mac_bins, labels = mac_bins_labels)
        return(d)
    }))
    counts <- aggregate(switches ~ wes_variant + mac_bin, data = dt, FUN = sum)
    tested <- aggregate(switches ~ wes_variant + mac_bin, data = dt, FUN = length)
    counts <- data.table(counts, tested = tested$switches)
    return(counts)
}

# get errors by gene (how many genes are perfectly phased)
aggregate_errors_by_gene <- function(files, variants){
    lst <- lapply(files, function(file){
        print(file)
        d <- fread_phased_sites(file)
        d$wes_variant <- d$locus %in% variants$locus
        d <- d[d$wes_variant == TRUE,]
        chrom <- gsub("chr","",unique(d$CHR))
        d <- get_switch_error_per_gene(pos = d$POS, switches = d$switches, chrom = chrom)
        d$chr <- paste0("chr",chrom)
        return(d)
    })
    return(lst)
}

# iterate over a matrix of switch errors and tested columns and append with
# binomial confidence intervals for the switch error rate
calc_binom_ci <- function(lst){
    counts <- do.call(rbind, lst)
    stopifnot(nrow(counts) > 0)
    counts_ci <- do.call(rbind, lapply(1:nrow(counts), function(i) Hmisc::binconf(counts$switches[i], counts$tested[i])))
    colnames(counts_ci) <- tolower(colnames(counts_ci))
    counts <- cbind(counts, counts_ci)
    return(counts)
}


main <- function(args){

    stopifnot(dir.exists(args$switch_dir))
    files <- list.files(args$switch_dir, pattern = args$switch_file)
    files <- file.path(args$switch_dir, files)
    print(files)
    stopifnot(length(files) > 0)
    stopifnot(length(files) > 21)  

    sites <- args$sites
    out_prefix <- args$out_prefix

    # read wes variants
    sites <- "/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list"
    variants <- fread(sites) 

    # get switch errors per gene
    lst_by_gene <- aggregate_errors_by_gene(files, variants)
    counts_by_gene <- do.call(rbind, lst_by_gene)
    counts_by_gene$wes_label <- "wes"
    
    out_gene <- paste0(out_prefix, "_gene.txt.gz")
    fwrite(counts_by_gene, out_gene, sep = "\t")
    write(paste("writing",out_gene), stderr())

    stop("stopped here")

    # create mac bin labels
    bins <- c(0,1,5,10,20,50,100,200,500,1000,2000,5000, 10000, Inf)
    labels <- unlist(lapply(2:length(bins), function(i){paste0(bins[i-1]+1,"-",bins[i])}))
    labels[labels == '1-1'] <- "singleton"

    # combine by chromosome
    lst_by_chrom <- aggregate_by_chrom(files, variants = variants)
    counts_by_chrom <- calc_binom_ci(lst_by_chrom)
    counts_by_chrom$wes_label <- ifelse(counts_by_chrom$wes_variant, "wes","array")

    out_chrom <- paste0(out_prefix, "_chrom.txt.gz")
    fwrite(counts_by_chrom, out_chrom, sep = "\t")
    write(paste("writing",out_chrom), stderr())

    # aggregate by mac bin and chromosome
    lst_by_mac_chrom <- aggregate_by_mac_chrom(files, bins, labels, variants)
    counts_by_mac_chrom <- calc_binom_ci(lst_by_mac_chrom)
    counts_by_mac_chrom$wes_label <- ifelse(counts_by_mac_chrom$wes_variant, "wes","array")

    out_mac_chrom <- paste0(out_prefix, "_mac_chrom.txt.gz")
    fwrite(counts_by_chrom, out_mac_chrom, sep = "\t")
    write(paste("writing",out_mac_chrom), stderr())

    # aggregate by mac bin
    dt_by_mac <- aggregate_by_mac(files, bins, labels, variants)
    counts_by_mac <- calc_binom_ci(list(dt_by_mac))
    counts_by_mac$wes_label <- ifelse(counts_by_mac$wes_variant, "wes","array")

    out_mac <- paste0(out_prefix, "_mac.txt.gz")
    fwrite(counts_by_mac, out_mac, sep = "\t")
    write(paste("writing",out_mac), stderr())


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--switch_dir", default=NULL, required = TRUE, help = "The directory containing chunks (to be searched recursively)")
parser$add_argument("--switch_file", default=NULL, required = FALSE, help = "Perform a subset (regex) based on a string")
parser$add_argument("--sites", default=NULL, required = TRUE, help = "path to quality controlled variant sites")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

