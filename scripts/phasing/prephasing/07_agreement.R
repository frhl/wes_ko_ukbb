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
    seed <- as.numeric(args$seed)
    sites <- args$sites
    input_path <- args$input_path    
    output_path <- args$output_path
    stopifnot(file.exists(sites))
    variants <- fread(sites)

    # get phased fil
    d <- fread(args$input_path)
    keep <- d$locus %in% variants$locus
    d <- d[keep,] 

    # keep phased sets with at least one rare variant
    bool_keep <- (!is.na(d$PP) & d$AC < 1001)
    ps_to_keep <- unique(d$PS_rb[bool_keep])
    dt <- d[d$PS_rb %in% ps_to_keep,]

    # use same bins as in S5 paper
    bins <- c(0,1,5,10,20,50,100,200,500,1000, Inf)
    labels <- unlist(lapply(2:length(bins), function(i){paste0(bins[i-1]+1,"-",bins[i])}))
    labels[labels == '1-1'] <- "singleton"
    dt$bin <- cut(dt$AC, breaks = bins, labels = labels)

    # sample 1000 random
    set.seed(seed)
    my_samples <- sample(unique(dt$s), n_samples, replace = FALSE)
    dt <- dt[dt$s %in% my_samples,]

    # eval 
    pp_50 <- eval_phase(eval_gt_agreement_by_bin(dt, labels[1:8]), labels)
    pp_60 <- eval_phase(eval_gt_agreement_by_bin(dt, labels[1:8]), labels)
    pp_70 <- eval_phase(eval_gt_agreement_by_bin(dt, labels[1:8]), labels)
    pp_75 <- eval_phase(eval_gt_agreement_by_bin(dt, labels[1:8]), labels)
    pp_80 <- eval_phase(eval_gt_agreement_by_bin(dt, labels[1:8]), labels)
    pp_85 <- eval_phase(eval_gt_agreement_by_bin(dt, labels[1:8]), labels)
    pp_90 <- eval_phase(eval_gt_agreement_by_bin(dt, labels[1:8]), labels)
    pp_95 <- eval_phase(eval_gt_agreement_by_bin(dt, labels[1:8]), labels)
    pp_99 <- eval_phase(eval_gt_agreement_by_bin(dt, labels[1:8]), labels)
    
    # annotate
    pp_50$pp <- 0.50
    pp_60$pp <- 0.60
    pp_70$pp <- 0.70
    pp_75$pp <- 0.75
    pp_80$pp <- 0.80
    pp_85$pp <- 0.85
    pp_90$pp <- 0.90
    pp_95$pp <- 0.95
    pp_99$pp <- 0.99

    # combine
    final <- rbind(pp_50, pp_60, pp_70, pp_75, pp_80, pp_85, pp_90, pp_95, pp_99)
    final$pp <- factor(final$pp)
    final$bin <- factor(final$bin, levels = labels) 
    fwrite(final, output_path, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, required = TRUE, help = "path to file with PS_rb, GT_rb and GT sites.")
parser$add_argument("--sites", default=NULL, required = TRUE, help = "path to quality controlled variant sites")
parser$add_argument("--output_path", default=NULL, required = TRUE, help = "")
parser$add_argument("--n_samples", default=10, required = TRUE, help = "")
parser$add_argument("--seed", default=52, required = TRUE, help = "")
args <- parser$parse_args()

main(args)

