#!/usr/bin/env Rscript


devtools::load_all("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/gwastools")
library(argparse)
library(data.table)

main <- function(args){

  print(args)

  stopifnot(file.exists(args$path_markers))
  stopifnot(file.exists(args$path_ac_by_phenotypes))
  stopifnot(file.exists(args$path_hash_by_phenotypes))

  min_mac <- as.numeric(args$min_mac)
  d_phenos <- fread(args$path_ac_by_phenotypes)
  d_marker <- fread(args$path_markers)

  # check integrity of phenotype
  phenotype <- args$phenotype
  stopifnot(phenotype %in% colnames(d_phenos))
  stopifnot("id" %in% colnames(d_phenos))
  stopifnot(min_mac >= 0)
  cols <- c("id", phenotype)

  ### filter by allele count ###

  # perform subset
  d_phenos <- d_phenos[,cols, with = FALSE]
  d_phenos <- d_phenos[ d_phenos[[phenotype]] >=  min_mac,]
  if (nrow(d_phenos) == 0) stop(paste0("No variants left after filtering on min_mac>=",min_mac))

  # check integrity of marker file
  stopifnot("consequence_category" %in% colnames(d_marker))
  stopifnot("rsid" %in% colnames(d_marker))

  # check how many annotations are present
  annotation <- unlist(strsplit(args$annotation, split = ','))
  annotation_found <- unlist(lapply(annotation, function(x) x %in% d_marker$consequence_category))
  if (sum(annotation_found) == 0) stop(paste("annotation", args$annotation, "was/were not in the data?"))
  if (sum(annotation_found) != length(annotation_found)) warning("Some csqs categories were not found!")

  # subset markers
  d_marker <- d_marker[d_marker$consequence_category %in% annotation,]
  if (nrow(d_marker) == 0) stop(paste0("No variants left after filtering by ", args$annotation))

  # combine the two files into a set of final variants
  markers <- intersect(
      d_marker$rsid, 
      d_phenos$id
  )

  if (length(markers) == 0) stop(paste("No variants left after combining", args$path_ac_by_phenotypes, "and", args$path_markers))

  ### filter by perfect LD ###

  # Prune markers in perfect LD
  d_hash <- fread(args$path_hash_by_phenotypes)
  stopifnot("id" %in% colnames(d_hash))
  d_hash <- d_hash[,cols, with = FALSE]
  d_hash <- d_hash[d_hash$id %in% markers, ]
  in_perfect_ld <- duplicated(d_hash[[phenotype]])
  d_hash <- d_hash[!in_perfect_ld, ]

  # count how many markers were affected
  n_perfect_ld <- sum(in_perfect_ld)
  if (n_perfect_ld > 0) write(paste("Pruned out", n_perfect_ld,"markers in perfect LD"), stdout())

  # ensure that only overlapping markers are kept
  markers <- intersect(
        markers,
        d_hash$id
  )

  if (length(markers) == 0) stop(paste("No variants left after combining", args$path_hash_by_phenotypes, "and", args$path_markers, "(Perfect LD pruning)"))
  
  # ensure variants follow input formatting required by saige
  markers <- gwastools::order_markers(markers)
  fwrite(data.table(x=markers), args$outfile, row.names = FALSE, col.names = FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "phenotype, single string.")
parser$add_argument("--path_markers", default=NULL, required = TRUE, help = "list of currently used variants seperated by comma")
parser$add_argument("--path_ac_by_phenotypes", default=NULL, required = TRUE, help = "file to allele count by phenotypes")
parser$add_argument("--path_hash_by_phenotypes", default=NULL, required = TRUE, help = "file containing dosage hashes for removing markers in perfect LD")
parser$add_argument("--min_mac", default=1, required = FALSE, help = "Allele count threshold, greater than or equal '>='")
parser$add_argument("--outfile", default=NULL, required = TRUE, help = "where should the subsetted markeres be written")
parser$add_argument("--annotation", default=NULL, required = TRUE, help = "Annotations to perform subset on")
args <- parser$parse_args()

main(args)









