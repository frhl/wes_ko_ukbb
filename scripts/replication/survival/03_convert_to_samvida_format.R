library(data.table)
library(argparse)

# Define the main function to perform the transformation and writing of files
main <- function(args){
  
  # Read the input file using fread from data.table
  data <- fread(args$input, header = TRUE, data.table = FALSE)
  
  # Iterate over the columns for each gene, starting from the second column (skipping 'eid')
  for (i in 2:ncol(data)) {
    # Create a data frame for the current gene
    gene_data <- data.frame(
      eid = data[[1]], # 'eid' is the first column
      chromosome = 'NA', # Placeholder for chromosome information
      ensembl_gene_id = names(data)[i], # The name of the current column is the gene ID
      hgnc_symbol = 'NA', # Placeholder for HGNC symbol
      annotation = data[[i]], # The data from the current column
      variants = 'NA' # Placeholder for variant information
    )
    
    # Define the output file name
    out_file <- sprintf("%s%s.tsv", args$out_prefix, names(data)[i])
    
    # Write the data frame to a file using fwrite from data.table
    fwrite(as.data.table(gene_data), out_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  }
}

# Set up argparse to handle command-line arguments
parser <- ArgumentParser()
parser$add_argument("--input", required=TRUE, help="Path to the input file")
parser$add_argument("--out_prefix", required=TRUE, help="Output file prefix for the resulting files")

# Parse the arguments
args <- parser$parse_args()

# Call the main function with the parsed arguments
main(args)

