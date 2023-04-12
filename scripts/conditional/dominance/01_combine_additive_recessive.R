devtools::load_all("utils/modules/R/gwastools")
library(argparse)
library(data.table)

# read VCF file
fread_vcf <- function(path){
    stopifnot(file.exists(path))
    command <- paste('zcat ', path, '| grep -v "##" ')
    vcf <- fread(cmd=command, sep = '\t')
    cols <- colnames(vcf)
    final_id_col <- suppressWarnings(max(which(is.na(as.numeric(colnames(vcf))))))
    metadata <- vcf[,1:final_id_col]
    genotypes <- vcf[,(final_id_col+1):ncol(vcf)]
    return(list(metadata=metadata, genotypes=genotypes))
}

# standardize genotypes
genotypes_to_integer <- function(genotypes){
    genotypes[genotypes=="."] <- NA
    genotypes[genotypes=="0.00000"] <- 0
    genotypes[genotypes=="1.00000"] <- 1
    genotypes[genotypes=="2.00000"] <- 2
    genotypes <- sapply(genotypes, as.numeric) # matrix
    return(genotypes)
}

# make header of VCF file
make_vcf_dosage_header <- function(chrom){
    vcf_format <- '##fileformat=VCFv4.2'
    vcf_entry <-  '##FORMAT=<ID=DS,Number=1,Type=Float,Description="">'
    vcf_filter <- '##FILTER=<ID=PASS,Description="All filters passed">"'
    vcf_i1 <- '##INFO=<ID=AC,Number=A,Type=Integer,Description="Knockout count multiplied by two">'
    vcf_i2 <- '##INFO=<ID=HASH,Number=A,Type=String,Description="Hash function applied to dosages">'
    vcf_contig <- paste0('##contig=<ID=',chrom,',length=81195210>')
    vcf_out <- paste(vcf_format, vcf_entry, vcf_filter, vcf_i1, vcf_i2, vcf_contig, sep = '\n')
    return(vcf_out)
}

main <- function(args){

    add <- fread_vcf(args$additive_path)
    rec <- fread_vcf(args$recessive_path)

    add$genotypes <- genotypes_to_integer(add$genotypes)
    rec$genotypes <- genotypes_to_integer(rec$genotypes)

    max_pos_rec <- max(rec$metadata$POS)
    add$metadata$POS <- add$metadata$POS + max_pos_rec

    # combine all metadata
    metadata <- rbind(rec$metadata, add$metadata)
    genotypes <- rbind(rec$genotypes, rec$genotypes)
    metadata_and_genotypes <- cbind(metadata, genotypes)

    # write list of additive variants
    add_variants <- paste(add$metadata$`#CHROM`, 
                          add$metadata$POS,
                          add$metadata$REF,
                          add$metadata$ALT,
                          sep= ":"
                          )

    out_table <- paste0(args$out_prefix, ".additive.txt")
    fwrite(data.table(x=add_variants), out_table, col.names=FALSE)

    # write header and append information
    vcf_out = make_vcf_dosage_header(args$chrom)
    outfile = paste0(args$out_prefix, ".vcf")
    writeLines(text = vcf_out, outfile)
    fwrite(genotypes, outfile, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)


}

parser <- ArgumentParser()
parser$add_argument("--additive_path", default=NULL, help = "")
parser$add_argument("--recessive_path", default=NULL, help = "")
parser$add_argument("--chrom", default=NULL, help = "")
parser$add_argument("--out_prefix", default=NULL, help = "")
args <- parser$parse_args()

main(args)


