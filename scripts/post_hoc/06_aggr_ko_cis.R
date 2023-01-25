
library(data.table)
library(argparse)

# mapping to HGNC symbol
bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensembl_to_hgnc <- bridge$hgnc_symbol
names(ensembl_to_hgnc) <- bridge$ensembl_gene_id


main <- function(args){
   
    # load all knockouts
    knockouts_paths <- list.files(args$knockout_dir, pattern = args$knockout_pattern, full.names = TRUE)    
    stopifnot(length(knockouts_paths)>0)

    # load knockouts
    lst_knockouts <- lapply(knockouts_paths, function(path){
        d <- fread(path)
        d$s <- as.character(d$s)
        d <- d[,c("gene_id","s","knockout","pKO")]
        d$chromosome <- stringr::str_extract(path, "chr[0-9]+")
        d$hgnc_symbol <- ensembl_to_hgnc[d$gene_id]
        return(d)
    })

    # names to keep
    cis <- "Compound heterozygote (cis)"
    chet <- "Compound heterozygote"
    homs <- "Homozygote"
    het <- "Heterozygote"

    # check for common knockouts
    dt_all <- do.call(rbind, lst_knockouts)
    counts <- data.table(table(dt_all$gene_id[dt_all$pKO > 0]))
    counts <- counts[rev(order(counts$N))]
    common_plofs <- counts[counts$N > 10000,]
    dt_all <- dt_all[!(dt_all$gene_id %in% common_plofs$V1),]
    print(paste("excluded",nrow(common_plofs),"common knockouts."))

    # create file with them 
    dt_out <- dt_all[(dt_all$knockout %in% c(chet, homs, cis)),]
    dt_out$is_chet <- dt_out$knockout %in% chet
    dt_out$is_hom <- dt_out$knockout %in% homs
    dt_out$is_cis <- dt_out$knockout %in% cis
   
    # write file 
    outfile <- paste0(args$out_prefix,".txt.gz")
    write(paste("writing", outfile), stderr())
    fwrite(dt_out, outfile)


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--knockout_dir", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--knockout_pattern", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


