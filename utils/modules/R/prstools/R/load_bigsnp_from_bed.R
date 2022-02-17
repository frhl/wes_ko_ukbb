#' @title load bigSNP from bed file
#' @param bed file path to plink file WITH .bed extension
#' @param verbose print misc stats
#' @export

load_bigsnp_from_bed <- function(bed, verbose = TRUE){
    
    # Read from bed/bim/fam, it generates .bk and .rds files.
    if (!file.exists.ext(bed, '.bk')) snp_readBed(bed)
    basename <- tools::file_path_sans_ext(bed)
    rds <- paste0(basename,'.rds')
    big_snp <- snp_attach(rds)

    # extract the SNP information from the genotype
    map <- big_snp$map[-3]
    names(map) <- c("chr", "rsid", "pos", "a1", "a0")

    # remove non-autosomes (e.g. chr8_KI270821v1_alt)
    autosomes <- paste0('chr',1:22)
    big_snp$map <- big_snp$map[big_snp$map$chr %in% autosomes,]
    G <- big_snp$genotypes

    # Rename the data structures
    CHR <- as.numeric(gsub('chr','',map$chr))
    POS <- map$pos

    # get the CM information from hapmap SNPs
    POS2 <- snp_as_genetic_position(CHR, POS, mapdir = "data/prs/1000-genomes-genetic-maps",genetic_map = 'hapmap')
    if (verbose){
        matches <- sum(POS2==0)/length(POS2) # hapmap has many less missing variants than omni
        write(paste0(100*(1-round(matches,5)),'% of variants are in genetic map (hapmap).'),stdout())
    }
    return(invisible(list(G = G, POS = POS, POS2 = POS2, map = map, fam = big_snp$fam)))
}




