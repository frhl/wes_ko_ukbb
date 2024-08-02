#!/usr/bin/env Rscript


devtools::load_all("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/gwastools")
library(argparse)
library(data.table)


filter_to_candidates <- function(dt, pheno, chr, csqs = c("pLoF", "damaging_missense"), verbose = TRUE){
    
    #' @param dt a data.frame with the columns, rsid, phenotype and consequence_category
    #' @param pheno phenotype to filter on
    #' @param chr chromosome to filter on
    #' @param csqs consequences to filter on (vector)
    
    # check input
    autosomes <- paste0("chr",1:22)
    stopifnot(chr %in% autosomes)
    stopifnot("rsid" %in% colnames(dt))
    stopifnot("phenotype" %in% colnames(dt))
    stopifnot("consequence_category" %in% colnames(dt))
    stopifnot(pheno %in% dt$phenotype)
    
    # build chromosome column
    dt$chr <- unlist(
        lapply(dt$rsid, function(x) unlist(strsplit(x, split = ':'))[1])
    )
    
    stopifnot(dt$chr %in% autosomes)
    n0 <- nrow(dt)
    
    # perform subsetting
    dt <- dt[dt$chr %in% chr,]
    dt <- dt[dt$phenotype %in% phenotype,]
    dt <- dt[dt$consequence_category %in% csqs,]
    
    # spit out some summaries
    n_markers <- length(unique(dt$rsid))
    n_genes <- length(unique(dt$ensembl_gene_id))
    msg <- paste("Note: Subsetting to", n_markers, "(unique) markers near", n_genes, "relevant gene(s) on", chr)
    if (verbose) write(msg, stderr())
        
    # change column names
    colnames(dt)[colnames(dt)=="rsid"] <- "id"
        
    return(dt)
    
}

annotate_by_allele_count <- function(dt, phenotype, min_mac){
    
    #' @param dt a data.table with columns, id, phenotypes where
    # each phenotype cell has the allele count for the phenotype
    #' @param phenotype phenotype column that is in dt
    #' @param min_mac integer above zero 
    
    # check input  
    stopifnot(phenotype %in% colnames(dt))
    stopifnot("id" %in% colnames(dt))
    stopifnot(min_mac >= 0)
    cols <- c("id", phenotype)
    
    
    # perform subset
    dt <- dt[,cols, with = FALSE]
    dt$keep_variant_ac <- dt[[phenotype]] >=  min_mac
    #dt <- dt[ d,]
    
    # change colnames
    colnames(dt)[1:2] <- c("id", "ac")
    dt$min_mac_filter <- min_mac
    
    return(dt)
    
}


annotate_by_perfect_ld <- function(dt, phenotype){
    
    #' @param dt a data.table with columns, id, phenotypes where
    # each phenotype cell is the hash of the dosage string
    #' @param phenotype phenotype column that is in dt
    
    # Prune markers in perfect LD
    stopifnot(phenotype %in% colnames(dt))
    stopifnot("id" %in% colnames(dt))
    cols <- c("id", phenotype)
    
    # perform subsets
    dt <- dt[,cols, with = FALSE]
    colnames(dt) <- c("id","hash")
    
    # find markers in perfect LD
    hash_duplicated <- duplicated(dt$hash)
    hash_in_perfect_ld <- dt$hash[hash_duplicated ]
    id_in_perfect_ld <- dt$id[dt$hash %in% hash_in_perfect_ld]
    n <- sum(hash_duplicated)
    
    if (n>0){
        
        # make mapping to markers in perfect LD
        lst <- lapply(hash_in_perfect_ld, function(hash){
            paste(dt$id[ dt$hash %in% hash], collapse = ';')
        })

        print(lst)

        hash_to_markers <- unlist(lst)
        names(hash_to_markers) <- hash_in_perfect_ld

        # filter to markers that are not in perfect LD
        dt$keep_variant_ld <- !hash_duplicated

        # keep track of the markers in perfect LD
        dt$perfect_ld_marker <- hash_to_markers[dt$hash]
        dt$perfect_ld_marker <- unlist(
            lapply(1:nrow(dt), function(idx){
                regex <- paste0(";?",dt$id[idx],";?")
                gsub(regex, "", dt$perfect_ld_marker[idx])
            })
        )
    } else {
        
        dt$keep_variant_ld <- TRUE
        dt$perfect_ld_marker <- NA
        
    }
    
    return(dt)  
}


main <- function(args){

  print(args)

  if (!file.exists(args$path_ac_by_phenotypes)) stop(paste(args$path_ac_by_phenotypes, "does not exist!"))
  if (!file.exists(args$path_hash_by_phenotypes)) stop(paste(args$path_hash_by_phenotypes, "does not exist!"))
  if (!file.exists(args$path_markers_by_gene)) stop(paste(args$path_markers_by_gene, "does not exist!"))

  chromosome <- args$chromosome
  phenotype <- args$phenotype
  min_mac <- as.numeric(args$min_mac)

  d_gene <- fread(args$path_markers_by_gene)
  d_phenos <- fread(args$path_ac_by_phenotypes)
  d_hash <- fread(args$path_hash_by_phenotypes)

  # here we need original phenotype when 'primary_care' is used.
  d_gene_filter <- filter_to_candidates(d_gene, phenotype, args$chromosome)
  markers_in_gene <- d_gene_filter$id

  # 
  if (!args$chromosome %in% d_gene_filter$chr) write(paste("Note: No significant genes for",phenotype,"on",chromosome), stdout())

  # filter hash and phenotype markers based on markers in gene
  d_phenos <- d_phenos[d_phenos$id %in% markers_in_gene,]
  d_hash <- d_hash[d_hash$id %in% markers_in_gene,]
  
  # annotate the data by various metrics
  d_ac_filter <- annotate_by_allele_count(d_phenos, phenotype, min_mac = 4)
  d_ld_filter <- annotate_by_perfect_ld(d_hash, phenotype)

  # merge the two filters
  mrg <- merge(d_ac_filter, d_ld_filter, by = 'id')
  stopifnot(nrow(mrg) == nrow(d_ac_filter))
  stopifnot(nrow(mrg) == nrow(d_ld_filter))

  # merge with gene ID and consequence
 
  ## Note, that the rsids may be duplicated if the same variant affects a two overlapping genes, e.g:
  # chr19:40778389:A:G    ENSG00000167578 damaging_missense   COPD_combined_primary_care  chr19
  # chr19:40778389:A:G    ENSG00000171570 damaging_missense   COPD_combined_primary_care  chr19 
  mrg <- merge(mrg, d_gene_filter, by = 'id', all.x = TRUE)
  mrg$chr <- NULL
  #mrg$phenotype <- NULL

  # markers to keep
  mrg$keep <- mrg$keep_variant_ac & mrg$keep_variant_ld
  
  if (sum(mrg$keep) > 0) {
    
    # create file containing the gene-specific variants
    outfile = paste0(args$outfile,".full")
    fwrite(mrg, outfile, sep = '\t', row.names = FALSE)

    # get final set of markers
    mrg <- mrg[mrg$keep,]
    markers <- mrg$id[!duplicated(mrg$id)]

    # create file containing the gene-specific variants
    outfile = paste0(args$outfile,".keep")
    fwrite(mrg, outfile, sep = '\t', row.names = FALSE)

    # ensure variants follow input formatting required by saige
    markers <- gwastools::order_markers(markers)
    fwrite(data.table(x=markers), args$outfile, row.names = FALSE, col.names = FALSE)

    # get counts for verbose message
    n_markers_start <- length(unique(d_gene_filter$id))
    n_markers_end <-  length(unique(mrg$id))
    n_genes <- length(unique(d_gene_filter$ensembl_gene_id))
    n_ac_ok <- sum(d_ac_filter$keep_variant_ac)
    n_ld_ok <- sum(!d_ld_filter$keep_variant_ld)
    n_dup_mark <- ifelse(length(unique(mrg$id)) != length(mrg$id),"Yes","No")
    n_ac_above_zero <- sum(d_ac_filter$ac > 0)
    n_ac_singletons <- sum(d_ac_filter$ac == 1)

    msg <- paste0(
      "# Phenotype: ",phenotype, "\n",
      "# Significant gene(s): ", n_genes, "\n",
      "# Maximum hail markers in gene(s): ", n_markers_start, "\n",
      "# Same marker in different gene(s): ", n_dup_mark, "\n",
      "# Markers in perfect LD: ", n_ld_ok, "\n",
      "# Markers with AC>=0: ", n_ac_above_zero, "\n",
      "# Markers with AC==1: ", n_ac_singletons, "\n",
      "# Markers with AC>=", min_mac,": ", n_ac_ok, "\n",
      "# Markers (post filtering of AC>=",min_mac," and LD) count: ", n_markers_end, ""
    )

    write(msg, stderr())

  } else {
    write("Note: No markers present after subsetting", stderr())
  }

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chromosome", default=NULL, required = TRUE, help = "Chromosome GRCh38, e.g. chr12")
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "phenotype, single string.")
parser$add_argument("--path_markers_by_gene", default=NULL, required = TRUE, help = "file to all markers within a given gene")
parser$add_argument("--path_ac_by_phenotypes", default=NULL, required = TRUE, help = "file to allele count by phenotypes")
parser$add_argument("--path_hash_by_phenotypes", default=NULL, required = TRUE, help = "file containing dosage hashes for removing markers in perfect LD")
parser$add_argument("--path_markers_in_chrom", default=NULL, required = FALSE, help = "file containing the markers in the chromosome for the current VCF")
parser$add_argument("--min_mac", default=1, required = FALSE, help = "Allele count threshold, greater than or equal '>='")
parser$add_argument("--outfile", default=NULL, required = TRUE, help = "where should the subsetted markeres be written")
parser$add_argument("--annotation", default=NULL, required = TRUE, help = "Annotations to perform subset on")
args <- parser$parse_args()

main(args)









