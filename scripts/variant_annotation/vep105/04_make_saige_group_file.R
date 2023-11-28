#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# Function to omit null values from a list
omit_null <- function(lst) {
    lst[!vapply(lst, is.null, logical(1))]
}

# Main function
main <- function(args) {
    # Read data from input file
    data <- fread(args$input_path)

    # any duplicated?
    n_in <- nrow(data)
    data <- data[!duplicated(data),]
    n_diff <- n_in - nrow(data)
    if (n_diff>0){warning(paste(n_diff,"duplicate records removed!"))}

    # Select and rename relevant columns
    selected_data <- data[, .(gene_id, varid, brava_csqs)]
    setnames(selected_data, c('gene', 'variant', 'annotation'))

    # Get unique genes
    genes <- unique(selected_data$gene)

    # Process each gene
    output <- lapply(genes, function(gene) {
        gene_variants <- selected_data$variant[selected_data$gene == gene]
        gene_annotations <- selected_data$annotation[selected_data$gene == gene]

        # Filter out NA values
        valid_indices <- !(is.na(gene_variants) | is.na(gene_annotations))
        gene_variants <- gene_variants[valid_indices]
        gene_annotations <- gene_annotations[valid_indices]

        # Define accepted annotation categories
        accepted_categories <- c('pLoF', 'damaging_missense', "other_missense", "synonymous", "non_coding")
        accepted <- gene_annotations %in% accepted_categories

        if (sum(accepted) > 0) {
            gene_variants <- gene_variants[accepted]
            gene_annotations <- gene_annotations[accepted]

            # Merge annotations if specified
            if (!is.null(args$merge_into_single_annotation)) {
                annotations_to_merge <- unlist(strsplit(args$merge_into_single_annotation, split = ","))
                combined_annotation <- paste0(annotations_to_merge, collapse = "_")
                stopifnot(all(annotations_to_merge %in% accepted_categories))
                gene_annotations[gene_annotations %in% annotations_to_merge] <- combined_annotation
            }

            # Create output rows
            variant_row <- paste(c(gene, 'var', gene_variants), collapse = " ")
            annotation_row <- paste(c(gene, 'anno', gene_annotations), collapse = " ")
            return(paste0(c(variant_row, annotation_row), collapse = '\n'))
        }
    })

    # Write output to file
    output <- omit_null(output)
    writeLines(paste(output, collapse = '\n'), args$output_path)
}

# Setup command line argument parser
parser <- ArgumentParser()
parser$add_argument("--merge_into_single_annotation", default = NULL, help = "Annotations to merge into a single category")
parser$add_argument("--input_path", default = NULL, help = "Path to the input file")
parser$add_argument("--output_path", default = NULL, help = "Path for the output file")
parser$add_argument("--delimiter", default = NULL, help = "Delimiter used in the input file")
args <- parser$parse_args()

# Run the main function
main(args)

