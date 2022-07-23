#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

  # read and check integratity of file
  stopifnot(file.exists(args$spa_file))
  d <- fread(args$spa_file, header = TRUE)
  cols <- colnames(d)
  stopifnot(nrow(d) > 0)
  
  # always use conditional P-valu column if present
  if (args$col_pvalue_cond %in% cols) {
    col_pvalue <- args$col_pvalue_cond
  } else if (args$col_pvalue %in% cols) {
    col_pvalue <- args$col_pvalue
  } else {
    stop(paste(args$col_pvalue,"/",args$col_pvalue_cond, "not in",args$spa_file))
  }
  
  # only keep the follow columns in the output file
  keep_cols <- c("CHR", "POS", "MarkerID", "Allele1", "Allele2",
                 "AC_allele2", "AF_Allele2", "BETA", "SE", "Tstat",
                 "var", "p.value", "Tstat_cond", "p.value_cond", "varT_vond",
                 "BETA_cond", "SE_cond") 
  
  # subset columns
  d <- d[,colnames(d) %in% keep_cols, with = FALSE]
  d$MARKER <- apply( d[ ,c("CHR", "POS", "Allele1", "Allele2") ] , 1 , paste , collapse = ":" )
  d$PVAL <- d[[col_pvalue]]
 
  # perform rows based on markers that pass thresholds
  p_cutoff <- as.numeric(args$p_cutoff)
  bool <- d$PVAL < p_cutoff
  d <- d[bool, ]
    
  # sort by most significant at BOTTOM (tail)
  d <- d[rev(order(d$PVAL)), ]

  #write("Example head and tail (n=2):", stdout())
  #print(head(d, n = 2))
  #print(tail(d, n = 2))
  # write file
  fwrite(d, args$out_file, sep = "\t", col.names = FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--spa_file", default=NULL, help = "?")
parser$add_argument("--out_file", default=NULL, help = "?")
parser$add_argument("--p_cutoff", default=NULL, help = "?")
parser$add_argument("--col_pvalue", default="p.value", help = "?")
parser$add_argument("--col_pvalue_cond", default="p.value_cond", help = "?")
args <- parser$parse_args()

main(args)

