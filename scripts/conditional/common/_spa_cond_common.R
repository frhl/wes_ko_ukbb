#!/usr/bin/env Rscript

devtools::load_all("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/gwastools")
library(argparse)
library(data.table)

main <- function(args){


    stopifnot(file.exists(args$infile))
    d <- fread(args$infile, header = FALSE)
    n <- nrow(d)
    
    # check that input are indexes 
    stopifnot(is.numeric(args$marker_col))
    if (!is.null(args$pheno_col)) stopifnot(grepl("[0-9]+",args$pheno_col))
    if (!is.null(args$annotation_col)) stopifnot(grepl("[0-9]+", args$marker_col))
    
    # ensure that they are actually integers/numerics
    args$pheno_col <- as.numeric(args$pheno_col)
    args$marker_col <- as.numeric(args$marker_col)
    args$annotation_col <- as.numeric(args$annotation_col)

    # subset by phenptype column
    if (!is.null(args$phenotype)){
        phenotypes <- d[[args$pheno_col]]
        print(phenotypes)
        if (!is.character(phenotypes)) stop(paste("Phenotype column", args$pheno_col, "does not contain characters!"))
        stopifnot(grepl('_',phenotypes))
        bool_pheno <- phenotypes %in% args$phenotype
        n_pheno <- n - sum(bool_pheno)
        d <- d[bool_pheno, ]

    }
    
    # subset by annotation column
    if (!is.null(args$annotation)){
        annotations <- d[[args$annotation_col]]
        bool_annotations <- annotations %in% args$phenotype
        stopifnot(is.character(annotations))
        stopifnot(grepl('_', annotations))
        n_annotation <- n - sum(bool_annotations)
        d <- d[bool_annotations, ]
    }
  
    
    # check markers format
    markers <- d[[args$marker_col]]
    regex <- "chr[0-9]+\\:[0-9]+\\:[a-zA-Z]+\\:[a-zA-Z]"
    if (any(!grepl(regex, markers))) stop("Markers are not formatted correctly!")
    if (length(markers) > 0) {
        markers <- gwastools::order_markers(markers)
        if (!is.null(args$outfile)){
            dout <- data.table(x=markers)
            outfile <- args$outfile
            write(paste0("Writing common markers to ",outfile), stderr())
            fwrite(dout, outfile, sep = '\t', col.names = FALSE)    
        } else {
            out <- paste(markers, collapse = ',')
            write(out, stdout())
        }
    } else {
        write("Note: No markers left after subsetting.", stderr())
    }


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--infile", default=NULL, help = "Table of markers")
parser$add_argument("--phenotype", default=NULL, help = "What phenotype subset should be made")
parser$add_argument("--annotation", default=NULL, help = "What annotation subset should be made")
parser$add_argument("--marker_col", default=2, help = "Index for column of markers")
parser$add_argument("--pheno_col", default=NULL, help = "Index for column of phenotypes")
parser$add_argument("--annotation_col", default=NULL, help = "Index for column for annotation")
parser$add_argument("--outfile", default=NULL, help = "Should a file be written? Default is write to stdout")
args <- parser$parse_args()

main(args)

