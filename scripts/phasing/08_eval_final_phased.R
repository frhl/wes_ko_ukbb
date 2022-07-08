#!/usr/bin/env Rscript

devtools::load_all('utils/modules/R/phasingtools')
library(argparse)
library(RColorBrewer)
library(stringr)

fread_phased_sites <- function(file, ...){
    
    # get details about chunks
    bname <- basename(file)
    chunk_current <- as.numeric(gsub("of","",stringr::str_extract(bname, "[0-9]+of")))
    chunk_final <- as.numeric(gsub("of","",stringr::str_extract(bname, "of[0-9]+")))
    method <- unlist(strsplit(bname, split = '_'))[1]
    phasing_region_size <- as.numeric(gsub("_prs","",stringr::str_extract(bname, "_prs[0-9]+")))
    phasing_overlap_size <- as.numeric(gsub("_pro","",stringr::str_extract(bname, "_pro[0-9]+")))
    max_phasing_region_size <- as.numeric(gsub("_mprs","",stringr::str_extract(bname, "_mprs[0-9]+")))
    
    # append to data.table
    d <- fread(file, ...)
    d$locus <- paste0(d$CHR,":",d$POS)
    d$chunk_current <- chunk_current
    d$chunk_final <- chunk_final
    d$method <- method
    d$phasing_region_size <- phasing_region_size
    d$phasing_overlap_size <- phasing_overlap_size
    d$max_phasing_region_size <- max_phasing_region_size
    return(d)
    
}

# aggregate switch errors by chromosome and minor allele frequency bins
aggregate_by_chrom <- function(files){
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

# same as above, but also aggregate my selected MAF bin
aggregate_by_chrom_and_maf_bin <- function(files, maf_bin){
    lst <- lapply(files, function(file){
        d <- fread_phased_sites(file)
        d$wes_variant <- d$locus %in% variants$locus
        d$maf_bin <- cut(d$MAF, breaks = cuts)
        counts <- aggregate(switches ~ wes_variant + maf_bin + CHR, data = d, FUN = sum)
        tested <- aggregate(switches ~ wes_variant + maf_bin + CHR, data = d, FUN = length)
        counts <- data.table(counts, tested = tested$switches)
        return(counts)
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

    # parser
    print(args)
    stopifnot(dir.exists(args$ligated_dir))
    stopifnot(dir.exists(dirname(args$out_prefix)))
    maf_bins <- c(1, 10^-(1:6))

    files <- list.files(args$ligated_dir, pattern = ".txt", full.names = TRUE)
    stopifnot(length(files) > 0)
    autosomes <- paste0("chr",1:22)
 
    lst_by_chrom <- aggregate_by_chromsome(files)
    counts_by_chrom <- calc_binom_ci(lst_by_chrom)

    # *** plot switch errors by chromosome ***
    pd <- position_dodge(1)
    p1 <- ggplot(counts_by_chrom,
           aes(
               y=factor(CHR, levels = autosomes),
               x=100*pointest,
               xmax = 100*upper,
               xmin = 100*lower,
               fill = factor(wes_variant)
           )) +
        geom_bar(stat = 'identity', position = pd, size = 1) +
        geom_errorbar(stat='identity', position = pd,width = 0.75) +
        geom_point(position = pd) +
        labs(fill = "WES variant") +
        xlab('Switch Errors (%)') + ylab('') +
        theme_bw()

    # *** plot switch errors by chrom and MAF ****
    lst_by_chrom_by_maf = aggregate_by_chrom_and_maf_bin(files, maf_bins) 
    counts_by_chrom_by_maf <- calc_binom_ci(lst_by_chrom_by_maf)
    counts_by_chrom_by_maf <- counts_by_chrom_by_maf[counts_by_chrom_by_maf$wes_variant == TRUE,]

    pd <- position_dodge(0.7)
    plt <- ggplot(counts_by_chrom_by_maf,
           aes(
               x=100*pointest,
               xmax = 100*upper,
               xmin = 100*lower,
               y = maf_bin,
               fill = factor(wes_variant)
           )) +
        #geom_bar(stat = 'identity', position = pd, size = 1) +
        geom_point(stat='identity', position = pd, size = 2) +
        geom_errorbar(stat='identity', position = pd,width = 0.75) +
        labs(fill = "WES variant") +
        ylab('Switch Errors (%)') + xlab('') +
        theme_bw() +
        facet_wrap(~factor(CHR, levels = autosomes)) +
        coord_cartesian(xlim=c(0, 2))

    # *** counts by maf bin only across chroms ***
    aggr_d <- count_by_chrom_by_maf
    aggr_switches <- aggregate(switches ~ wes_variant + maf_bin, data = aggr_d, FUN = sum)
    aggr_tested <- aggregate(tested ~ wes_variant + maf_bin, data = aggr_d, FUN = sum)
    aggr_counts <- data.table(aggr_switches, tested = aggr_tested$tested)
    aggr_counts <- calc_binom_ci(aggr_counts)

    pd <- position_dodge(1)
    plt <- ggplot(aggr_counts,
           aes(
               x=100*pointest,
               xmax = 100*upper,
               xmin = 100*lower,
               y = maf_bin,
               fill = factor(wes_variant)
           )) +
        geom_bar(stat = 'identity', position = pd) +
        geom_point(stat='identity', position = pd, size = 2) +
        geom_errorbar(stat='identity', position = pd,width = 0.50) +
        labs(fill = "WES variant") +
        xlab('Switch Errors (%)') + ylab('MAF bin') +
        theme_bw()

    # need to

    #write(paste0("writing to",out_p1), stdout())
    #ggsave(p1, out_p1, width = args$img_width, height = args$img_height) 

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--ligated_dir", default=NULL, required = TRUE, help = "The directory containing chunks (to be searched recursively)")
parser$add_argument("--sites", default=NULL, required = TRUE, help = "path to quality controlled variant sites")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--img_width", default=8, help = "Where should the results be written?")
parser$add_argument("--img_height", default=6, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

