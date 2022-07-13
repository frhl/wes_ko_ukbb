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
    if (!is.null(args$files_regex)) files <- files[grepl(args$files_regex, files)]
    stopifnot(length(files) > 0)
    autosomes <- paste0("chr",1:22)
    print(files)
    

    lst_by_chrom <- aggregate_by_chromsome(files)
    counts_by_chrom <- calc_binom_ci(lst_by_chrom)
    counts_by_chrom$wes_label <- ifelse(counts_by_chrom$wes_variant, "Whole Exome Sequencing","Genotyping Array")

    # *** plot switch errors by chromosome ***
    p1 <- ggplot(counts_by_chrom,
       aes(
           y=factor(CHR, levels = autosomes),
           x=100*pointest,
           xmax = 100*upper,
           xmin = 100*lower,
           color = factor(wes_label)
       )) +
    theme_bw() +
    geom_pointrange() + 
    labs(color = "") +
    xlab('% Switch Error Rate (95% CI)') + ylab('') +
    scale_color_d3('category20c', limits=NULL) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    theme(
        legend.position = "top",
        axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"),
        axis.title.x = element_text(margin=ggplot2::margin(t=10)),
        axis.title.y = element_text(margin=ggplot2::margin(r=10)),
        plot.title = element_text(hjust=0.5)
    ) 

    out_p1 <- paste0(args$out_prefix, "_ser_by_chrom")
    out_p1_img <- paste0(out_p1, ".png")
    out_p1_txt <- paste0(out_p1, ".txt.gz")
    ggsave(p1, out_p1_img, width = 10, height = 8)
    fwrite(counts_by_chrom, out_p1_txt, sep = "\t")

    # *** plot switch errors by chrom and MAF ****
    lst_by_chrom_by_maf = aggregate_by_chrom_and_maf_bin(files, maf_bins) 
    counts_by_chrom_by_maf <- calc_binom_ci(lst_by_chrom_by_maf)
    counts_by_chrom_by_maf <- counts_by_chrom_by_maf[counts_by_chrom_by_maf$wes_variant == TRUE,]
    counts_by_chrom_by_maf$wes_label <- ifelse(counts_by_chrom_by_maf$wes_variant, "Whole Exome Sequencing","Genotyping Array")
    
    p2 <- ggplot(counts_by_chrom_by_maf,
           aes(
               y=maf_bin,
               x=100*pointest,
               xmax = 100*upper,
               xmin = 100*lower,
               color = factor(wes_label)
           )) +
        theme_bw() +
        geom_pointrange() + 
        labs(color = "") +
        xlab('% Switch Error Rate (95% CI)') + ylab('') +
        scale_color_d3('category20c', limits=NULL) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
        facet_wrap(~factor(CHR, levels = autosomes)) +
        coord_cartesian(xlim=c(0, 2)) +
        theme(
            legend.position = "top",
            axis.text=element_text(size=10),
            axis.title=element_text(size=10,face="bold"),
            axis.title.x = element_text(margin=ggplot2::margin(t=10)),
            axis.title.y = element_text(margin=ggplot2::margin(r=10)),
            plot.title = element_text(hjust=0.5)
        ) 

    out_p2 <- paste0(args$out_prefix, "_ser_by_maf_chrom")
    out_p2_img <- paste0(out_p2, ".png")
    out_p2_txt <- paste0(out_p2, ".txt.gz")
    ggsave(p2, out_p2_img, width = 10, height = 8)
    fwrite(counts_by_chrom, out_p2_txt, sep = "\t")


    # *** counts by maf bin across all chromosomes ***
    aggr_d <- count_by_chrom_by_maf
    aggr_switches <- aggregate(switches ~ wes_variant + maf_bin, data = aggr_d, FUN = sum)
    aggr_tested <- aggregate(tested ~ wes_variant + maf_bin, data = aggr_d, FUN = sum)
    aggr_counts <- data.table(aggr_switches, tested = aggr_tested$tested)
    aggr_counts <- calc_binom_ci(aggr_counts)
    aggr_counts$wes_label <- ifelse(aggr_counts$wes_variant, "Whole Exome Sequencing","Genotyping Array")

    p3 <- ggplot(aggr_counts,
           aes(
               y=maf_bin,
               x=100*pointest,
               xmax = 100*upper,
               xmin = 100*lower,
               color = factor(wes_label)
           )) +
        theme_bw() +
        geom_pointrange() + 
        labs(color = "") +
        xlab('% Switch Error Rate (95% CI)') + ylab('') +
        scale_color_d3('category20c', limits=NULL) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
        theme(
            legend.position = "top",
            axis.text=element_text(size=10),
            axis.title=element_text(size=10,face="bold"),
            axis.title.x = element_text(margin=ggplot2::margin(t=10)),
            axis.title.y = element_text(margin=ggplot2::margin(r=10)),
            plot.title = element_text(hjust=0.5)
        ) 

    out_p3 <- paste0(args$out_prefix, "_ser_by_maf")
    out_p3_img <- paste0(out_p1, ".png")
    out_p3_txt <- paste0(out_p1, ".txt.gz")
    ggsave(p3, out_p3_img, width = 8, height = 6)
    fwrite(counts_by_chrom, out_p3_txt, sep = "\t")


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--ligated_dir", default=NULL, required = TRUE, help = "The directory containing chunks (to be searched recursively)")
parser$add_argument("--sites", default=NULL, required = TRUE, help = "path to quality controlled variant sites")
parser$add_argument("--files_regex", default=NULL, required = FALSE, help = "Perform a subset (regex) based on a string")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

