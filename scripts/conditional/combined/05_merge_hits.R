
library(argparse)
library(data.table)


main <- function(args){

    # first, get markers that we have tested by conditionion on everything
    sig_genes <- fread(args$path_sig_genes)
    traits_to_test <- fread(args$path_header, header = FALSE)$V1
    sig_genes <- sig_genes[sig_genes$trait %in% traits_to_test,]

    # combine list of significant genes 
    d <- rbindlist(lapply(1:nrow(sig_genes), function(idx){
        phenotype <- sig_genes$trait[idx]
        gene <- sig_genes$gene[idx]
        regex <- paste0(phenotype,".+",gene,"$")
        file_path <- list.files(args$target_dir, pattern = regex, full.names = TRUE)
        prs <- grepl("locoprs", file_path)
        if (length(file_path) == 0){
            write(paste0("path '",file_path,"' (",gene,"/",phenotype,") could not be found. Skipping."), stderr())
        } else {
            # open file and subset to the gene
            d <- fread(file_path)
            d <- d[d$MarkerID %in% gene,]
            d$phenotype <- phenotype
            d$prs <- prs
            d$filepath <- file_path
            return(d)
        }
    }))

    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(d, outfile, sep = "\t", quote=FALSE)

}

parser <- ArgumentParser()
parser$add_argument("--target_dir", default=NULL, help = "")
parser$add_argument("--path_header", default=NULL, help = "")
parser$add_argument("--path_sig_genes", default=NULL, help = "")
parser$add_argument("--out_prefix", default=NULL, help = "")
args <- parser$parse_args()

main(args)


