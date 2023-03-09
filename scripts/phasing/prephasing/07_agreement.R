#!/usr/bin/env Rscript

devtools::load_all('utils/modules/R/phasingtools')
library(argparse)
library(data.table)
library(Hmisc)

bin_conf <- function(counts){
    counts_ci <- do.call(rbind, lapply(1:nrow(counts), function(i) Hmisc::binconf(counts$match[i], counts$total[i])))
    colnames(counts_ci) <- tolower(colnames(counts_ci))
    counts <- cbind(counts, counts_ci)                                   
    return(counts)
}

eval_phase <- function(full, labels){
    sums <- aggregate(match ~ bin, data = full, FUN = sum)
    lens <- aggregate(match ~ bin, data = full, FUN = length)
    colnames(lens)[2] <- "total"
    mrg <- merge(sums, lens)
    mrg <- mrg[match(labels, mrg$bin), ]
    mrg <- mrg[!is.na(mrg$bin),]
    mrg$errors <- mrg$total - mrg$match
    mrg <- bin_conf(mrg)
    return(mrg)
}
                                       

main <- function(args){

    # get WES sites
    n_samples <- as.numeric(args$n_samples)
    pp_cutoff <- as.numeric(args$pp_cutoff)
    seed <- as.numeric(args$seed)
    sites <- args$sites
    mac_bin <- args$mac_bin
    input_path <- args$input_path    
    out_prefix <- args$out_prefix
    stopifnot(file.exists(sites))
    variants <- fread(sites)

    # get phased fil
    d <- fread(args$input_path)
    stopifnot("MAC" %in% colnames(d))
    stopifnot("PP" %in% colnames(d))
    keep <- d$locus %in% variants$locus
    d <- d[keep,] 
    d$MAC <- as.numeric(d$MAC)

    # bins
    bins <- c(0,1,5,10,20,50,100,200,500,1000,2000,5000,10000, Inf)
    
    # keep phased sets with at least one rare variant
    bool_keep <- (!is.na(d$PP) & d$MAC < (bins[length(bins)-1]+1))
    ps_to_keep <- unique(d$PS_rb[bool_keep])
    dt <- d[d$PS_rb %in% ps_to_keep,]
    stopifnot(nrow(dt) > 0)

    # use same bins as in S5 paper
    labels <- unlist(lapply(2:length(bins), function(i){paste0(bins[i-1]+1,"-",bins[i])}))
    labels[labels == '1-1'] <- "singleton"
    dt$bin <- cut(dt$MAC, breaks = bins, labels = labels)
   
    # ensure that selected mac_bin is available.
    stopifnot(mac_bin %in% labels)
    labels_to_run <- mac_bin

    # sample 1000 random
    set.seed(seed)
    my_samples <- sample(unique(dt$s), n_samples, replace = FALSE)
    dt <- dt[dt$s %in% my_samples,]

    # NA's are due to well phased genotypes
    dt$PP[is.na(dt$PP)] <- 1

    # eval gt agreement
    gt_agreement <- phasingtools::eval_gt_agreement_by_bin(dt, mac_bin, pp_cutoff)
    print(head(gt_agreement))
    
    # return genotype agreement
    the_bin <- gsub("-", "_", mac_bin)
    outfile <- paste0(out_prefix,".",the_bin,".txt.gz")
    fwrite(gt_agreement, outfile, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, required = TRUE, help = "path to file with PS_rb, GT_rb and GT sites.")
parser$add_argument("--sites", default=NULL, required = TRUE, help = "path to quality controlled variant sites")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "out prefix for files.")
parser$add_argument("--n_samples", default=10, required = TRUE, help = "how many samples should be sampled?")
parser$add_argument("--seed", default=52, required = TRUE, help = "seed for selecting random samples")
parser$add_argument("--pp_cutoff", default=NULL, required = TRUE, help = "What phasing confidence threshold should be cut of?")
parser$add_argument("--mac_bin", default=NULL, required = TRUE, help = "Should a particular mac_bin be evaluated? Faster than running all.")
args <- parser$parse_args()

main(args)

