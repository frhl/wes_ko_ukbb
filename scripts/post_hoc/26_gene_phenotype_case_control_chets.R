
devtools::load_all("utils/modules/R/gwastools")
source("scripts/post_hoc/utils.R")
library(data.table)
library(argparse)


main <- function(args){

    # get things to keep
    co_occurence_file = args$co_occurence_file
    out_prefix = args$out_prefix

    # go over all chroms
    lst <- list()
    for (chrom in 1:22){
        path <- gsub("CHR",chrom, co_occurence_file)
        co <- fread(path)
        lst[[chrom]] <- co
    }
    d <- rbindlist(lst)

    # aggregate by variant ID
    aggr1 <- setDT(aggregate(chet~g+phenotype+type, data=d, FUN=sum))
    aggr2 <- setDT(aggregate(hom~g+phenotype+type, data=d, FUN=sum))
    aggr3 <- setDT(aggregate(cis~g+phenotype+type, data=d, FUN=sum))
    keys <- c("g","phenotype", "type")
    setkeyv(aggr1, keys)
    setkeyv(aggr2, keys)
    setkeyv(aggr3, keys)
    
    #  combine across all gene-traits
    mrg <- merge(merge(aggr1, aggr2), aggr3)
    mrg$gene_trait <- paste0(mrg$phenotype,":",mrg$g)

    outfile <- paste0(out_prefix,".txt.gz")
    fwrite(mrg, outfile, sep = "\t")
 
}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--co_occurence_file", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


