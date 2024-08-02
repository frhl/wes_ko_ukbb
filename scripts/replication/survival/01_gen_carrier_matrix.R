#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){

    input_file <- args$input_file
    sample_file <- args$sample_file
    out_prefix <- args$out_prefix
    n_min_ko <- as.numeric(args$n_min_ko)
    genes_per_block <- args$genes_per_block
    set_genes_to_test <- args$set_genes_to_test
    chrom <- args$chrom

    d <- fread(input_file, header=FALSE)
    print(head(d))
    colnames(d) <- c("s", "contig", "gene_id", "config", "ds", "variant")
    d$variant <- NULL
    d$ds <- NULL
   
    # get samples 
    samples <- fread(sample_file, header=FALSE)$V1
    d <- d[d$s %in% samples,]
    if (!is.null(chrom)) d <- d[d$contig %in% chrom,]

    # get counts
    gene_freq <- data.table(table(d$gene_id, d$config))
    gene_freq <- gene_freq[gene_freq$N != 0,]

    # convert to wide
    gene_freq <- dcast(V1~V2, data = gene_freq)
    colnames(gene_freq)[1] <- "gene_id"
    gene_freq[is.na(gene_freq)] <- 0
    gene_freq$ko <- gene_freq$chet + gene_freq$hom
    gene_freq <- gene_freq[rev(order(gene_freq$ko)),]

    # We test if there is at least 5 KOs
    if (!is.null(set_genes_to_test)){
        stopifnot(file.exists(set_genes_to_test))
        genes_to_test <- fread(set_genes_to_test, header=FALSE)$V1
        print("Here are genes to test before subset:")
        print(genes_to_test)
        print("We are checking against this dataset (that 'genes to test' are in here):")
        print(head(gene_freq))
        genes_to_test <- unique(genes_to_test[genes_to_test %in% gene_freq$gene_id])
        if (length(genes_to_test) == 0) stop("No genes left to test in data!")
        print(paste("Manually gene-testing selected:", paste0(genes_to_test, collapse=",")))
    } else {
        gene_freq <- gene_freq[gene_freq$ko >= n_min_ko]
        genes_to_test <- gene_freq$gene_id
    }

    print(paste("Combining", length(genes_to_test), "genes with at least one chet or hom."))
    if ( length(genes_to_test) > 0) {
        lst <- lapply(genes_to_test, function(g){
            print(g)
            # remove redundant info
            d_by_g <- d[d$gene_id %in% g,]
            non_wt <- length(d_by_g$s)
            d_by_g$contig <- NULL
            d_by_g$gene_id <- NULL
            # append missing samples (assume these are WT) 
            missing_samples <- samples[!samples %in% d_by_g$s]
            d_by_g <- rbind(d_by_g, data.table(s=missing_samples, config="wt"))
            d_by_g <- d_by_g[order(d_by_g$s),]
            colnames(d_by_g)[2] <- g
            return(d_by_g)    
        })
        combined <- Reduce(merge, lst)

        # write out all genes as a single block
        if (is.null(genes_per_block)){
            out <- paste0(out_prefix, ".txt.gz")
            fwrite(combined, out, sep="\t")
        # write out genes in seperate blocks
        } else {
            genes <- colnames(combined)[-1]
            n_genes <- length(genes)
            shift <- as.numeric(genes_per_block)
            total_blocks <- ceiling(n_genes / shift)
            block_file <- paste0(out_prefix, ".blocks.txt")
            col_offset <- 1
            for (idx in 1:total_blocks){
                  # get slices
                  start_idx <- (idx-1)*shift+1+col_offset
                  end_idx <- idx*shift+col_offset
                  end_idx <- min(end_idx, ncol(combined))
                  indicies <- start_idx:end_idx
                  # subset columns
                  cols <- c('s', colnames(combined[,..indicies]))
                  block <- combined[,..cols]
                  id <- paste0("b", idx, "of", total_blocks)
                  # write block
                  out <- paste0(out_prefix,".",id,".txt.gz")
                  cat(out, file=block_file, sep="\n", append=TRUE)
                  fwrite(block, out, sep="\t")
            } 
        }
    } else {
        stop("No genes left to test!")
    } 
}

parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = FALSE, help = "Chromosome")
parser$add_argument("--set_genes_to_test", default=NULL, required = FALSE, help = "Chromosome")
parser$add_argument("--genes_per_block", default=NULL, required = FALSE, help = "Number of genes per block")
parser$add_argument("--n_min_ko", default=1, required = FALSE, help = "Minimum number of KOs")
parser$add_argument("--input_file", default=NULL, required = TRUE, help = "Input_file")
parser$add_argument("--sample_file", default=NULL, required = TRUE, help = "sample_file")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()


main(args)

