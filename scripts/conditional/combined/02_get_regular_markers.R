#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(digest)

fread_vcf <- function(path){
    stopifnot(file.exists(path))
    command <- paste('zcat ', path, '| grep -v "##"')
    vcf <- fread(cmd=command, sep = '\t')
    cols <- colnames(vcf)
    final_id_col <- suppressWarnings(max(which(is.na(as.numeric(colnames(vcf))))))
    metadata <- vcf[,1:final_id_col]
    genotypes <- vcf[,(final_id_col+1):ncol(vcf)]
    return(list(metadata=metadata, genotypes=genotypes))
}


main <- function(args){

    # read mac/maf count (which we can't read from VEP file due to additional QC/filtering)
    m <- fread(args$path_markers)
    # read variant file with variant to gene
    v <- fread(args$path_worst_csq_by_gene_canonical)
    m <- m[m$MAC > 0,]
    m <- m[!duplicated(m),]
    v <- v[v$varid %in% m$rsid,] 
    v <- v[,c("locus","alleles","worst_csq_by_gene_canonical.gene_id")]
    v <- v[!duplicated(v),]
    # combine the two files
    mv <- merge(m, v, by = c("locus","alleles"))
    colnames(mv)[6] <- "ensembl_gene_id"
    genes <- unique(mv$ensembl_gene_id)

    # read vcf and meta-data
    l <- fread_vcf(args$path_vcf)

    # get genes and their corresponding variants
    pseudo_rsids <- lapply(genes, function(g){
        return(mv$rsid[mv$ensembl_gene_id %in% g])
    })

    genotypes <- l$genotypes
    metarow <- l$metadata
    row_gene <- toupper(metarow$ID)
    row_marker <- paste0(metarow$`#CHROM`,":",metarow$POS,":",metarow$REF,":",metarow$ALT)

    # ensure correct format
    genotypes[genotypes=="."] <- NA
    genotypes[genotypes=="0.00000"] <- 0
    genotypes[genotypes=="1.00000"] <- 1
    genotypes[genotypes=="2.00000"] <- 2
    genotypes <- sapply(genotypes, as.numeric) # matrix

    # get AC and hash
    pseudo_AC <- apply(genotypes, 1, function(gts) sum(gts, na.rm = TRUE))
    pseudo_hash <- apply(genotypes, 1, function(gts) digest(gts, algo="xxhash32"))

    # get overview
    final <- data.table(
        chrom = args$chrom,
        ensembl_gene_id = row_gene,
        psuedo_marker = row_marker,
        pseudo_AC = pseudo_AC,
        pseudo_hash = pseudo_hash
    )

    # get genes and their corresponding variants by category
    categories <- unique(mv$consequence_category)
    pseudo_csqs <- rbindlist(lapply(genes, function(g){
        ds <- data.table(table(mv$consequence_category[mv$ensembl_gene_id %in% g]))
        ds <- data.table(t(ds))
        colnames(ds) <- as.character(ds[1,])
        ds <- ds[2,]
        # ugly hack to add missing categories
        categories_missing <- categories[!categories %in% colnames(ds)]
        for (cat in categories_missing){ds[[cat]] <- NA}
        ds <- ds[,..categories]
        ds$ensembl_gene_id <- g
        return(ds)
    }))
    colnames(pseudo_csqs)[1:2] <- paste0("n_",colnames(pseudo_csqs)[1:2])
    
    # merge and write file
    final_with_csqs <- merge(final, pseudo_csqs,  all.x = TRUE)
    outfile = paste0(args$out_prefix, ".txt.gz")
    write(paste("Succes! Writing to", outfile), stdout())
    fwrite(final_with_csqs, outfile, quote = FALSE, sep = '\t', col.names = TRUE)
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_worst_csq_by_gene_canonical", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_markers", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_vcf", default=NULL, required = TRUE, help = "")
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "where should the subsetted markeres be written")
args <- parser$parse_args()

main(args)









