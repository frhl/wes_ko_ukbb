#' @title load bigSNP from bed file
#' @param bed file path to plink file WITH .bed extension
#' @param verbose print misc stats
#' @export

load_bigsnp_from_bed <- function(bed, verbose = TRUE){
    
    # Read from bed/bim/fam, it generates .bk and .rds files.
    #tmp <- tempfile(tmpdir = "data/tmp/tmp-data")
    #on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
    

    if (!file.exists.ext(bed, '.bk')) snp_readBed(bed)
    basename <- tools::file_path_sans_ext(bed)
    rds <- paste0(basename,'.rds')
    
    # if the bk file exists but not the rds, then re-read SNP
    if (!file.exists(rds)) {
       write("rds file not available. Retrying snp_readBED", stderr())
       bk <- paste0(basename, ".bk")
       unlink(bk)
       snp_readBed(bed)
    }
    
    # attach to the file to current session
    big_snp <- snp_attach(rds)

    # extract the SNP information from the genotype
    autosomes <- paste0('chr',1:22)
    map <- big_snp$map[-3]
    names(map) <- c("chr", "rsid", "pos", "a1", "a0")
    
    #  check for odd chromsome contigs (usually remnants of liftover)
    odd_contig <- unique(map$chr[!map$chr %in% autosomes])
    if (length(odd_contig) > 0){
        odd_contig <- paste0(odd_contig, collapse = '|')
        warning(paste0("non standard chromosome contig excluded: ", odd_contig))
    }

    # remove non-standard chromosome contigs
    map <- map[map$chr %in% autosomes,]
    big_snp$map <- big_snp$map[big_snp$map$chr %in% autosomes,]
    G <- big_snp$genotypes

    # Rename the data structures
    CHR <- as.numeric(gsub('chr','',map$chr))
    POS <- map$pos

    # get the CM information from hapmap SNPs
    POS2 <- snp_as_genetic_position(CHR, POS, mapdir = "data/prs/1000-genomes-genetic-maps",genetic_map = 'hapmap')
    if (verbose){
        matches <- sum(POS2==0)/length(POS2) # hapmap has many less missing variants than omni
        write(paste0(100*(1-round(matches,5)),'% of variants are in genetic map (hapmap).'),stderr())
    }
    return(invisible(list(G = G, POS = POS, POS2 = POS2, map = map, fam = big_snp$fam, bigsnp = big_snp)))
}




