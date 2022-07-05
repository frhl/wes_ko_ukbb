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

main <- function(args){

    # parser
    print(args)
    stopifnot(dir.exists(args$master_chunk_dir))
    stopifnot(dir.exists(dirname(args$out_prefix)))

    
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

    # plot chunk by switch error rate stratified by variant type (exome or genotyping)
    pd <- position_dodge(0.7)
    p1 <- ggplot(counts,
           aes(
               y=100*pointest,
               ymax = 100*upper,
               ymin = 100*lower,
               x = factor(chunk_current), #factor(CHR, levels = autosomes),
               fill = factor(wes_variant)
           )) +
        geom_bar(stat = 'identity', position = pd, size = 1) +
        geom_errorbar(stat='identity', position = pd,width = 0.75) +
        labs(fill = "WES variant") +
        ylab('Switch Errors (%)') + xlab('') +
        theme_bw() +
        facet_wrap(~factor(CHR, levels = autosomes))

    out_p1 <- paste0(args$out_prefix, "_chunk_by_ser.png")
    write(paste0("writing to",out_p1), stdout())
    ggsave(p1, out_p1, width = args$img_width, height = args$img_height) 

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

