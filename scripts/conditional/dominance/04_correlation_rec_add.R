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


main <- function(args){

    # read file and get markers 
    chrom <- args$chrom
    dom <- fread_vcf(args$input_path)
    dom$genotypes <- genotypes_to_integer(dom$genotypes)
    recessive_mates <- dom$metadata$ID[grepl(dom$metadata$ID,pattern="^ENSG")]
    
    # find mate and calculate correlation
    final_by_chrom <- rbindlist(lapply(recessive_mates, function(recessive_mate){
        additive_mate <- gsub("ENSG","addensg", recessive_mate)
        idx_recessive <- which(dom$metadata$ID == recessive_mate)
        idx_additive <- which(dom$metadata$ID == additive_mate)
        ds_recessive <- dom$genotypes[idx_recessive,]
        ds_additive <- dom$genotypes[idx_additive,]
        correlation <- cor(ds_recessive, ds_additive)
        out <- data.frame(recessive_mate, additive_mate, correlation, chrom)
        return(out)
    }))

    # write file
    out_file <- paste0(args$out_prefix, ".txt")
    write(paste("writing", out_file), stdout())
    fwrite(final_by_chrom, out_file, sep="\t")


}

parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "")
parser$add_argument("--chrom", default=NULL, help = "")
parser$add_argument("--out_prefix", default=NULL, help = "")
args <- parser$parse_args()

main(args)


