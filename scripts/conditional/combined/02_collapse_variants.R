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

# make header of VCF file
make_vcf_dosage_header <- function(chrom){
    vcf_format <- '##fileformat=VCFv4.2'
    vcf_entry <-  '##FORMAT=<ID=DS,Number=1,Type=Float,Description="">'
    vcf_filter <- '##FILTER=<ID=PASS,Description="All filters passed">"'
    vcf_contig <- paste0('##contig=<ID=',chrom,',length=81195210>')
    vcf_out <- paste(vcf_format, vcf_entry, vcf_filter, vcf_contig, sep = '\n')
    return(vcf_out)
}

# make vcf-like rows for dosage entries
make_vcf_dosage_rows <- function(chrom, positions, marker){
    return(data.table(
      "#CHROM" = chrom,
      POS = positions,
      ID = marker,
      REF = 'R',
      ALT = 'A',
      QUAL = '.',
      FILTER = '.',
      INFO = '.',
      FORMAT = 'DS'
    ))
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

    # get genes and their corresponding variants
    pseudo_rsids <- lapply(genes, function(g){
        return(mv$rsid[mv$ensembl_gene_id %in% g])
    })

    # read vcf and meta-data
    l <- fread_vcf(args$path_vcf)
    genotypes <- l$genotypes
    write(paste("Setting dosages..",args$out_prefix), stderr())

    # ensure correct format
    genotypes[genotypes=="."] <- NA
    genotypes[genotypes=="0.00000"] <- 0
    genotypes[genotypes=="1.00000"] <- 1
    genotypes[genotypes=="2.00000"] <- 2
    genotypes <- sapply(genotypes, as.numeric) # matrix

    write(paste("Making pseudo dosages..",args$out_prefix), stderr())
    # get pseuedo gts
    pseudo_gts <- lapply(pseudo_rsids, function(rsids){
        idx <- which(l$metadata$ID %in% rsids)
        g <- genotypes[idx,]
        print(nrow(g))
        print(ncol(g))
        pseudo_gt <- apply(g, 2, function(x) max(x, na.rm = TRUE))
        return(pseudo_gt)
    })

    # LD hash 
    write(paste("Making hashes..",args$out_prefix), stderr())
    pseudo_ac_hash <- lapply(pseudo_gts, function(gts){
        return(data.table(
            ac = sum(gts, na.rm = TRUE),
            hash = digest(gts, algo="xxhash32") 
        ))
    })

    # combine
    dt_pseudo_gts <- do.call(rbind, pseudo_gts)
    dt_ac_hash <- do.call(rbind, pseudo_ac_hash)

    # get pseudo meta data
    max_pos <- max(l$metadata$POS)
    n_genes <- length(genes)
    rows <- make_vcf_dosage_rows(args$chrom, max_pos+(1:n_genes), paste0(genes,"_rare"))
    final <- cbind(rows, dt_pseudo_gts)

    # (1) write header of VCF
    vcf_out = make_vcf_dosage_header(args$chrom)
    outfile = paste0(args$out_prefix, ".vcf")
    writeLines(text = vcf_out, outfile)

    # (2) append with permuted data
    fwrite(final, outfile, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)

    # (3) write marker IDs
    genes <- genes
    markers <- paste(rows$`#CHROM`, rows$POS, rows$REF, rows$ALT, sep = ':')
    cond_markers <- cbind(genes, markers, dt_ac_hash)
    outmarkers = paste0(args$out_prefix, ".markers")
    fwrite(final, outmarkers, quote = FALSE, sep = '\t', col.names = TRUE)
    write(paste("Done!..",args$out_prefix), stderr())

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









