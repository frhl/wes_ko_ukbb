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
    stopifnot(dir.exists(args$master_chunk_dir))
    stopifnot(dir.exists(dirname(args$out_prefix)))
    stopifnot(file.exists(dirname(args$sites)))

    # we are restricting to chromosome 20-22 for this comparison
    files <- list.files(args$master_chunk_dir, pattern = ".txt", full.names = TRUE, recursive = TRUE)
    files <- files[grepl("chr2[0-2]", files)]
    variants <- fread(args$sites)
    
    if (!is.null(args$files_regex)) files <- files[grepl(args$files_regex, files)]
    stopifnot(length(files) > 0)
    autosomes <- paste0("chr",1:22)
    
    lst_by_method <- lapply(files, function(f){
        d <- data.table(do.call(rbind, lapply(f, fread_phased_sites)))
        d$wes_variant <- d$locus %in% variants$locus
        counts <- aggregate(switches ~ wes_variant + chunk_current + CHR, data = d, FUN = sum)
        tested <- aggregate(switches ~ wes_variant + chunk_current + CHR, data = d, FUN = length)
        counts <- data.table(counts, tested = tested$switches)
        counts$method <- ifelse(grepl("eagle", f), "Eagle2","SHAPEIT4")
        return(counts)
    })

    # get CIs for estimates
    counts_by_method <- calc_binom_ci(lst_by_method)
    pd <- position_dodge(0.7)
    p1 <- ggplot(counts_by_method,
           aes(
               y=switches,
               x = factor(chunk_current), #factor(CHR, levels = autosomes),
               fill = factor(wes_variant)
           )) +
        geom_bar(stat = 'identity', position = pd, color = 'black') +
        ylab('Switch Error Rate (%)') + xlab('') +
        labs(fill = "WES variant") +
        theme_bw() +
        facet_grid(factor(CHR, levels = autosomes)~method)
 
    out_p1 <- paste0(args$out_prefix,'_eagle_shapeit4.png')       
    out_d1 <- paste0(args$out_prefix,'_eagle_shapeit4.txt.gz')       
    write(paste0("writing to",out_p1), stdout())
    ggsave(p1, out_p1, width = args$img_width, height = args$img_height) 
    fwrite(counts_by_method, out_d1, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--master_chunk_dir", default=NULL, required = TRUE, help = "The directory containing chunks (to be searched recursively)")
parser$add_argument("--sites", default=NULL, required = TRUE, help = "path to quality controlled variant sites")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--img_width", default=5, help = "Where should the results be written?")
parser$add_argument("--img_height", default=8, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

