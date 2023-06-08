devtools::load_all("utils/modules/R/gwastools")
library(argparse)
library(data.table)

# load get current genes used
source("scripts/post_hoc/utils.R")

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

    # check if all mate-pairs are there
    recessive_mates <- dom$metadata$ID[grepl(dom$metadata$ID,pattern="^ENSG")]
    additive_mates <- dom$metadata$ID[grepl(dom$metadata$ID,pattern="^addensg")]

    # check missing mates
    missing_mates_idx <- which(!recessive_mates %in% gsub("addensg","ENSG",additive_mates))
    missing_mates <- recessive_mates[missing_mates_idx]

    # check if missing mates are in expected genes (remember we removed some pseudogenes)
    # but did not change the original recessive coding for these few genes.
    if (length(missing_mates) > 0){
        expected_genes <- read_ukb_wes_kos("pLoF_damaging_missense", chrom=chrom)
        expected <- sum(missing_mates %in% expected_genes$gene_id)
        if (sum(expected) > 0) stop(paste("Additive mate for gene",expected_genes,"was expected but not found! "))
    }

    # get subset that is used in recessive encoding
    mates_to_test <- gsub("addensg","ENSG",additive_mates)
    stopifnot(all(mates_to_test %in% recessive_mates))

    # find mate and calculate correlation
    final_by_chrom <- rbindlist(lapply(mates_to_test, function(recessive_mate){
        additive_mate <- gsub("ENSG","addensg", recessive_mate)
        idx_recessive <- which(dom$metadata$ID == recessive_mate)
        idx_additive <- which(dom$metadata$ID == additive_mate)
        #write(additive_mate, stderr())
        #write(recessive_mate, stderr())
        ds_recessive <- dom$genotypes[idx_recessive,]
        ds_additive <- dom$genotypes[idx_additive,]
        # ensure that there is always an additive mate
        if (length(idx_recessive) == 0) stop(paste("Missing the recessive mate", recessive_mate))
        if (length(idx_additive) == 0) stop(paste("Missing the additive mate for", recessive_mate))
        stopifnot(length(ds_recessive) == length(ds_additive))
        # check that recessive is truly a subset of additive and calulcate correlation
        ds_additive_wo_hets <- ds_additive
        ds_additive_wo_hets[ds_additive_wo_hets==1] <- 0
        correlation <- cor(ds_recessive, ds_additive)
        warn_nomatch <- as.integer(!all(ds_additive_wo_hets == ds_recessive))
        if (warn_nomatch) {
            warning(paste("Recessive encoding is not a subset of additve encoding:", args$chrom, recessive_mate))
            correlation <- NA               
        } 
        out <- data.frame(recessive_mate, additive_mate, correlation, chrom, warn_nomatch)    
    }))
    
    # counts
    final_by_chrom$mates_tested <- length(mates_to_test)
    final_by_chrom$recessive_mates <- length(recessive_mates)
    final_by_chrom$additive_mates <- length(additive_mates)

    # write file
    out_file <- paste0(args$out_prefix, ".txt")
    write(paste("writing", out_file), stdout())
    fwrite(final_by_chrom, out_file, sep="\t",na="NA")


}

parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "")
parser$add_argument("--chrom", default=NULL, help = "")
parser$add_argument("--out_prefix", default=NULL, help = "")
args <- parser$parse_args()

main(args)


