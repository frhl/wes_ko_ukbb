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

    # get marker to gene
    marker_to_gene <- mv$ensembl_gene_id
    names(marker_to_gene) <- mv$rsid
    marker_to_csqs <- mv$consequence_category
    names(marker_to_csqs) <- mv$rsid

    # read vcf and meta-data
    l <- fread_vcf(args$path_vcf)
    genotypes <- l$genotypes
    metarow <- l$metadata
    row_marker <- as.character(paste0(metarow$`#CHROM`,":",metarow$POS,":",metarow$REF,":",metarow$ALT))
    row_gene <- as.character(marker_to_gene[row_marker])
    row_csqs <- as.character(marker_to_csqs[row_marker])

    # ensure correct format
    genotypes[genotypes=="."] <- NA
    genotypes[genotypes=="0.00000"] <- 0
    genotypes[genotypes=="1.00000"] <- 1
    genotypes[genotypes=="2.00000"] <- 2
    genotypes <- sapply(genotypes, as.numeric) # matrix

    # get AC and hash
    pseudo_AC <- apply(genotypes, 1, function(gts) sum(gts, na.rm = TRUE))
    pseudo_hash <- apply(genotypes, 1, function(gts) digest(gts, algo="xxhash32"))

    # combine everything 
    final <- data.table(
        ensembl_gene_id = row_gene,
        marker = row_marker,
        AC = pseudo_AC,
        hash = pseudo_hash,
        csqs = row_csqs
    )

    # merge and write file
    outfile = paste0(args$out_prefix, ".txt.gz")
    write(paste("Succes! Writing to", outfile), stdout())
    fwrite(final, outfile, quote = FALSE, sep = '\t', col.names = TRUE)
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









