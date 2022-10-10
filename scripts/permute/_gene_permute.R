#!/usr/bin/env Rscript

library(argparse)
library(data.table)


# hash dosage as a string. This function simply concatenates all the dosgae counts
# into a long string, and subsequently applies a hash to this string. 
hash_dosage <- function(DS_matrix, algo = "xxhash32"){
    require(digest)
    dosage_string <- unlist(apply(DS_matrix, 1, function(x) as.character(paste(x, collapse = '-'))))
    the_hash <- unlist(lapply(dosage_string, function(x) digest(x, algo=algo)))
    return(the_hash)
}

# get probability of knockout given number of heterozygotes
naive_knockout_p <- function(hets){
    if (hets > 1){
        return(1 - (2*((1/2)^hets)))
    } else {
        return(0)
    }
}

# fast way to calculate probability of being KO
calc_knockout_p_fast <- function(d){
    
    # hom ref are always zero
    d$P <- 0
    # hom alt are always one
    d$P[d$hom_alt_n > 1] <- 1
    # depends on count for phased hets
    index <- which((d$phased_het > 1) & (d$hom_alt_n == 0))
    for (idx in index){
        hets <- d$phased_het[idx]
        d$P[idx] <- naive_knockout_p(hets)
    }
    
    return(d$P)
    
}

# simple method to shuffle knockouts
shuffle_knockouts <- function(d){
    
    # header of d
    # gene_id   s   unphased_het    phased_het  hom_alt_n   het pTKO
    # <chr> <int>   <int>   <int>   <int>   <int>   <dbl>
    # ENSG00000027644   1000028 0   0   0   0   0
    # ENSG00000027644   1000034 0   0   0   0   0
    # ENSG00000027644   1000087 0   0   0   0   0
    
    n <- nrow(d)
    p <- calc_knockout_p_fast(d)
    B <- rbinom(n = n, size = 1, prob = p)

    return(B)
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

# aggregate information from dosage matrix to create INFO rows
calc_info <- function(DS_matrix){
    stopifnot(ncol(DS_matrix) > 1)
    stopifnot(nrow(DS_matrix) > 1)
    # create INFO row by comining strings column wise
    d_info <- data.table(
        allele_count = rowSums(DS_matrix),
        hashes = hash_dosage(DS_matrix)    
    )
    # append identifiers as specified in header
    d_info$allele_count <- paste0("AC=",d_info$allele_count,";")
    d_info$hashes <- paste0("HASH=",d_info$hashes)
    info_rows <- apply(d_info, 1, paste, collapse = "")
    return(info_rows)
}


# make a random string of n characters and numbers
random_string <- function(n){
    return(paste0(sample(c(letters, LETTERS), n, replace = TRUE), collapse = ""))
} 

# make vcf-like rows for dosage entries
make_vcf_dosage_rows <- function(chrom, positions, marker, use_random_alleles = TRUE){
    rows <- length(positions)
    refs <- unlist(ifelse(use_random_alleles, list(replicate(rows, random_string(4))), "A0"))
    alts <- unlist(ifelse(use_random_alleles, list(replicate(rows, random_string(4))), "A1"))
    return(data.table(
      "#CHROM" = chrom, 
      POS = positions,
      ID = marker,
      REF = refs,
      ALT = alts,
      QUAL = '.',
      FILTER = '.',
      INFO = '.',
      FORMAT = 'DS'
    ))
}

# read real variant data in long format (exported from hail)
# note: need position_last argument, since make_tabix will otherwise complain
# that the positions are not sorted.
format_real_variant_long_to_wide <- function(dt, position_last = 20000){
    
    stopifnot("s" %in% colnames(dt))
    stopifnot("locus" %in% colnames(dt))
    stopifnot("rsid" %in% colnames(dt))
    stopifnot("alleles" %in% colnames(dt))
    
    mapping <- dt[,c("locus","alleles","rsid")]
    mapping <- mapping[!duplicated(mapping),]

    # create mapping rows that are to be combined with actual dosages
    alleles <- as.data.frame(do.call(rbind, strsplit(gsub('(")|(\\])|(\\[)','',mapping$alleles), split = ',')))
    colnames(alleles) <- c("REF","ALT")
    mapping <- cbind(mapping, alleles)
    mapping$positions_recoded <- (position_last+1):(position_last+nrow(mapping))
    mapping$positions <- gsub("chr[0-9]+\\:", "",mapping$locus)
    mapping$marker <- paste0(mapping$locus, ":", mapping$REF, ":", mapping$ALT)
    mapping$chroms <- stringr::str_extract(mapping$locus, "chr[0-9]+")
    stopifnot(length(unique(mapping$positions)) == length(mapping$positions)) # can't handle SNPs at same pos
    mapping_rows <- data.table(
        "#CHROM" = mapping$chroms,
          POS = mapping$positions,
          ID = mapping$rsid,
          REF = mapping$REF,
          ALT = mapping$ALT,
          QUAL = '.',
          FILTER = '.',
          INFO = '.',
          FORMAT = 'DS',
          locus = mapping$locus
    )
    
    # go from long to wide format
    dt <- dt[,c('locus','s','DS')]
    dt <- data.table::dcast(locus~s, data = dt, value.var = "DS")
    
    # match mapping rows with dt rows
    #new_index <- match(dt$locus, mapping$locus)
    new_index <- match(mapping$locus, dt$locus)
    mapping_rows$locus <- NULL
    #mapping_rows <- mapping_rows[new_index,]
    dt <- dt[new_index,]
    
    #return(cbind(mapping_rows, dt))
    return(list(rows = mapping_rows, dosages = dt))
    
}


main <- function(args){

    # print(args)

    autosomes <- paste0("chr",1:22)
    stopifnot(file.exists(args$input_path))
    stopifnot(!is.na(as.numeric(args$permutations)))
    stopifnot(!is.null(args$permutations))
    stopifnot(!is.null(args$vcf_id))
    stopifnot(args$chrom %in% autosomes)
    stopifnot(!is.null(args$enable_cond_pipeline))
    # conditional pipeline
    if (args$enable_cond_pipeline) {
      cpath = args$input_path_cond_genotypes  
      if (!file.exists(cpath)){
        stop(paste("Can't find conditional markers at: ", cpath))
      }
    }

    # seed for reproducibility
    seed <- as.numeric(args$seed)
    set.seed(seed)

    # replicate knockout
    n <- as.numeric(args$permutations)
    d <- fread(args$input_path)
    stopifnot(nrow(d) > 0)
    reps <- replicate(n, shuffle_knockouts(d))
    rownames(reps) <- d$s
    reps <- data.table(t(reps))

    # convert to dosage
    dosage <- reps * 2

    # only keep non-probabilistic knockouts
    if (args$only_non_prob_ko) dosage[dosage < 2] <- 0


    # load real conditioning variants, i.e. the actual
    # dosages/genotypes of the variants that we would like
    # to condition on
    if (args$enable_cond_pipeline) {
      cond_dt <- fread(args$input_path_cond_genotypes)
      cond_dt$chr <- stringr::str_extract(cond_dt$locus, "chr[0-9]+")
      cond_dt <- cond_dt[cond_dt$chr %in% args$chrom]
      n_real_markers <- length(unique(cond_dt$locus))
    } else {
      n_real_markers <- 0
    }

    # if there are conditiONIng markers available include them downstream.
    if (n_real_markers > 0){

        # how many markers were found?
        write(paste("Note:",n_real_markers, "real marker(s) found. These will be included as unshuffled in permuted VCF."),stdout())

        # ensure that samples are overlapping
        sample_overlap <- unique(intersect(cond_dt$s, d$s))

        # subset dosage matrix (with permuted phased)
        rows <- make_vcf_dosage_rows(args$chrom, 1:n, args$vcf_id)
        dosage <- dosage[,colnames(dosage) %in% sample_overlap, with = FALSE]

        # subset real dosage matrix (with actual calls/DS)
        cond_dt <- cond_dt[cond_dt$s %in% sample_overlap,]
        cond_lst <- format_real_variant_long_to_wide(cond_dt, n)

        # get long format
        cond_rows <- cond_lst$rows
        cond_dosage <- cond_lst$dosage

        # match columns
        cond_dosage$locus <- NULL
        cond_dosage <- cond_dosage[,colnames(dosage), with = FALSE]

        # combine columns and rows. Note: that rbind order
        # matters here when using tabix!
        stopifnot(ncol(rows) == ncol(cond_rows))
        combined_dosages <- rbind(dosage, cond_dosage) 
        combined_meta <- rbind(rows, cond_rows)
        final <- cbind(combined_meta, combined_dosages)
        
        # calculate some info stats
        sds <- unlist(apply(combined_dosages, 1, function(x) sd(x, na.rm = TRUE)))

        # calculate INFO column
        final$INFO <- calc_info(combined_dosages)

        # how many real markers have missing dosages
        n_real_count <- nrow(cond_dosage)
        n_real_miss <- sum(is.na(apply(cond_dosage, 1, sd)))
        write(paste0("Note: ", n_real_miss, " of ", n_real_count, " real markers have one or more missing dosages."), stdout())

        # debugging - are SNPs monoprhic and thus
        # the resulting matrix not invertible?
        #cond_dosage_sd <- apply(cond_dosage, 1, sd)
        #cond_dosage_af <- apply(cond_dosage, 1, mean)
        #cond_rows$sd <- cond_dosage_sd
        #cond_rows$af <- cond_dosage_af

    }  else {

        rows <- make_vcf_dosage_rows(args$chrom, 1:n, args$vcf_id)
        rows_dosage <- cbind(rows, dosage)
        final <- rows_dosage
        final$INFO <- calc_info(dosage)
        sds <- unlist(apply(dosage, 1, sd))
    }



    # Sometimes markers with zero AC are crated,
    # let's remove them before entering SAIGE.
    if (any(is.na(sds))) stop("Some standard devations are NA! Something went wrong with shuffle")
    if (args$remove_invariant_markers){
        bool_invariant <- sds == 0
        n_invariant <- sum(bool_invariant)
        print(n_invariant)
        if (n_invariant > 0){
            final <- final[!bool_invariant,]
            write(paste("[_gene_permute.R]: Removed", n_invariant, "invariant markers."), stderr())
        }
    }

    # (1) write header of VCF
    vcf_out = make_vcf_dosage_header(args$chrom)
    outfile = paste0(args$out_prefix, ".vcf")
    writeLines(text = vcf_out, outfile)

    # (2) append with permuted data
    fwrite(final, outfile, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, help = "chromosome")
parser$add_argument("--input_path", default=NULL, help = "path to the input")
parser$add_argument("--input_path_cond_genotypes", default=NULL, help = "path to the file of dosages/genotypes")
parser$add_argument("--permutations", default=NULL, help = "number of times the gene should be permuted")
parser$add_argument("--only_non_prob_ko", action="store_true", default=FALSE, help = "Only keep knockouts of phased heterozygotes.")
parser$add_argument("--remove_invariant_markers", action="store_true", default=FALSE, help = "Remove markers with AC == 0.")
parser$add_argument("--enable_cond_pipeline", action="store_true", default=FALSE, help = "Allow the use of conditional markers")
parser$add_argument("--seed", default=NULL, help = "seed for randomizer")
parser$add_argument("--vcf_id", default="GENE", help = "Substitute for rsid")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

