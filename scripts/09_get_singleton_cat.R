#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

    autosomes <- paste0(1:22)[20:21]
    lst <- list()
    singleton_lst <- lst()
    for (chr in autosomes){
        path_vep <- gsub("CHR", chr, args$vep_file)
        path_alt <- gsub("CHR", chr, args$in_file)
        stopifnot(file.exists(path_vep))
        stopifnot(file.exists(path_alt))
        
        # read files
        alt <- fread(path_alt, header=FALSE, sep=" ")
        colnames(alt) <- c("SAMPLE", "ID", "GT", "PP", "AC", "AN", "AF")
        vep <- fread(path_vep)

        print(head(alt))
        print(head(vep))

        # map signletons to csqs
        csqs_map <- vep$brava_csqs
        names(csqs_map) <- vep$varid
        singletons <- alt[alt$AC==1,]
        singletons$csqs <- csqs_map[singletons$ID]
    
        # save singletons
        singleton_lst[[chr]] <- data.table(ID=singletons$ID)

        # get matrix
        transposed_d <- t(data.table(table(singletons$csqs)))
        final <- data.table(transposed_d)
        colnames(final) <- transposed_d[1, ]
        final <- final[2, ]
        
        # convert category
        final$pLoF <- as.numeric(final$pLoF)
        final$damaging_missense <- as.numeric(final$damaging_missense)
        final$other_missense <- as.numeric(final$other_missense)
        final$synonymous <- as.numeric(final$synonymous)
        final$non_coding <- as.numeric(final$non_coding)

        # save to list
        final$chrom <- chr
        lst[[chr]] <- final
    }

    # combine final
    combine <- do.call(rbind, lst)
     
    # shwo the categories
    categories <- c("pLoF", "damaging_missense", "other_missense", "non_coding", "synonymous")
    for (cat in categories){
        count <- sum(combine[[cat]])
        print(paste(count, cat, "variants found across", length(unique(combine$chrom)), "chroms"))
    }
 
    # save main count
    out <- paste0(args$out_prefix, ".csqs.overview.txt.gz")
    write(paste("writing to", out), stdout()) 
    fwrite(combine, out, sep="\t") 

    # also save singletons
    singletons_combined <- do.call(rbind, singleton_lst) 
    out <- paste0(args$out_prefix, ".singletons.txt.gz")
    write(paste("writing to", out), stdout()) 
    fwrite(singletons_combined, out, sep="\t") 


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--vep_file", default=NULL, help = "")
parser$add_argument("--in_file", default=NULL, help = "")
parser$add_argument("--out", default=NULL, help = "")
args <- parser$parse_args()

main(args)

