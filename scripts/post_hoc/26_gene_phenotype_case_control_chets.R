
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
    for (chrom in 2:22){
        path <- gsub("CHR",chrom, co_occurence_file)
        co <- fread(path)
        lst[[chrom]] <- co
    }
    d <- rbindlist(lst)

    # aggregate by variant ID
    aggr1 <- setDT(aggregate(chet~g+phenotype+type, data=d, FUN=sum))
    aggr2 <- setDT(aggregate(hom~g+phenotype+type, data=d, FUN=sum))
    aggr3 <- setDT(aggregate(cis~g+phenotype+type, data=d, FUN=sum))
    aggr4 <- setDT(aggregate(het~g+phenotype+type, data=d, FUN=sum))
    keys <- c("g","phenotype", "type")
    setkeyv(aggr1, keys)
    setkeyv(aggr2, keys)
    setkeyv(aggr3, keys)
    setkeyv(aggr4, keys)
    
    #  combine across all gene-traits
    mrg <- merge(merge(merge(aggr1, aggr2), aggr3), aggr4)
    mrg$gene_trait <- paste0(mrg$phenotype,":",mrg$g)

    outfile <- paste0(out_prefix,".txt.gz")
    fwrite(mrg, outfile, sep = "\t")

    # convert to wide
    d <- mrg
    setkeyv(d, c("g","phenotype"))

    # count up cases
    d1 <- d[d$type == "cases",]
    d1 <- d1[,c("g","phenotype","chet","hom", "cis", "het")]
    colnames(d1)[colnames(d1) == "chet"] <- "N_ko_case.chetonly"
    colnames(d1)[colnames(d1) == "hom"] <- "N_ko_case.homonly"
    colnames(d1)[colnames(d1) == "cis"] <- "N_cis_case"
    colnames(d1)[colnames(d1) == "het"] <- "N_het_case"

    # count up cases and controls
    d2 <- d[d$type == "all",]
    d2 <- d2[,c("g","phenotype","chet","hom", "cis", "het")]
    colnames(d2)[colnames(d2) == "chet"] <- "N_ko.chetonly"
    colnames(d2)[colnames(d2) == "hom"] <- "N_ko.homonly"
    colnames(d2)[colnames(d2) == "cis"] <- "N_cis"
    colnames(d2)[colnames(d2) == "cis"] <- "N_het"
    # merge all
    mrg <- merge(d1, d2)
    outfile <- paste0(out_prefix,".wide.txt.gz")
    fwrite(mrg, outfile, sep = "\t")


}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--co_occurence_file", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


