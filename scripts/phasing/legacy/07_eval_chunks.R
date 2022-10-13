#!/usr/bin/env Rscript

devtools::load_all('utils/modules/R/phasingtools')
library(argparse)
library(RColorBrewer)
library(stringr)
library(ggsci)

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

main <- function(args){

    # parser
    print(args)
    stopifnot(dir.exists(args$master_chunk_dir))
    stopifnot(file.exists(args$sites))
    stopifnot(dir.exists(dirname(args$out_prefix)))
    autosomes <- paste0("chr",1:22)

    # load files   
    files <- list.files(args$master_chunk_dir, recursive = TRUE, pattern = ".txt", full.names = TRUE)
    files <- files[grepl("shapeit", files)]
    print("this is the header")
    print(head(files))
    variants <- fread(args$sites)

    lst <- lapply(autosomes, function(chr){
        
        # subset by chromosome (if we read in all at the same time it takes too long)
        files_chr <- files[grepl(paste0(chr,"-"), files)]
        d <- data.table(do.call(rbind, lapply(files_chr, fread_phased_sites)))
        d$wes_variant <- d$locus %in% variants$locus
        
        # get counts
        counts <- aggregate(switches ~ wes_variant + chunk_current + CHR, data = d, FUN = sum)
        tested <- aggregate(switches ~ wes_variant + chunk_current + CHR, data = d, FUN = length)
        counts <- data.table(counts, tested = tested$switches)
        return(counts)
    })  

    # setup counts for bar charts
    counts <- do.call(rbind, lst)
    counts_ci <- do.call(rbind, lapply(1:nrow(counts), function(i) Hmisc::binconf(counts$switches[i], counts$tested[i])))
    colnames(counts_ci) <- tolower(colnames(counts_ci))
    counts <- cbind(counts, counts_ci)
    counts$wes_label <- ifelse(counts$wes_variant, "Whole Exome Sequencing","Genotyping Array")

    # plot chunk by switch error rate stratified by variant type (exome or genotyping)i
    pd <- position_dodge(0.7)
    p1 <- ggplot(counts,
       aes(
           x=factor(chunk_current),
           y=100*pointest,
           ymax = 100*upper,
           ymin = 100*lower,
           fill = factor(wes_label)
       )) +
    theme_bw() +
    geom_bar(stat = 'identity', position = pd, size = 1) +
    geom_errorbar(stat='identity', position = pd,width = 0.75) +
    labs(color = "") +
    ylab('% Switch Error Rate (95% CI)') + xlab('Phasing chunks') +
    scale_color_d3('category20c', limits=NULL) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    facet_wrap(~factor(CHR, levels = autosomes)) +
    theme(
        legend.position = "top",
        axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"),
        axis.title.x = element_text(margin=ggplot2::margin(t=10)),
        axis.title.y = element_text(margin=ggplot2::margin(r=10)),
        plot.title = element_text(hjust=0.5)
    ) 

    h = as.numeric(args$img_height)
    w = as.numeric(args$img_width) 
    out_p1 <- paste0(args$out_prefix, "_chunks_by_ser.png")
    out_d1 <- paste0(args$out_prefix, "_chunks_by_ser.txt.gz")
    write(paste0("writing ",h,"x", w," png to",out_p1), stdout())
    fwrite(counts, out_d1, sep = "\t") 
    ggsave(out_p1, p1, width = w, height = h)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--master_chunk_dir", default=NULL, required = TRUE, help = "The directory containing chunks (to be searched recursively)")
parser$add_argument("--sites", default=NULL, required = TRUE, help = "path to quality controlled variant sites")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--img_width", default=8, help = "Where should the results be written?")
parser$add_argument("--img_height", default=6, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

