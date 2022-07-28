#!/usr/bin/env Rscript

library(argparse)
library(data.table)

z_to_log10p <- function(z) {
      abs((pnorm(-abs(z),log.p=TRUE)+log(2))/log(10))
}

main <- function(args){

  # read and check integratity of file
  stopifnot(file.exists(args$spa_file))
  d <- fread(args$spa_file, header = TRUE)
  cols <- colnames(d)
  stopifnot(nrow(d) > 0)
  
  # always use conditional P-valu column if present
  if (args$col_pvalue_cond %in% cols) {
    col_pvalue <- args$col_pvalue_cond
    col_tstat <- "Tstat_c"
    col_beta <- "BETA_c"
    col_se <- "SE_c"
  } else if (args$col_pvalue %in% cols) {
    col_pvalue <- args$col_pvalue
    col_tstat <- "Tstat"
    col_beta <- "BETA"
    col_se <- "SE"
  } else {
    stop(paste(args$col_pvalue,"/",args$col_pvalue_cond, "not in",args$spa_file))
  }
  
  # only keep the follow columns in the output file
  keep_cols <- c("CHR", "POS", "MarkerID", "Allele1", "Allele2",
                 "AC_allele2", "AF_Allele2", "BETA", "SE", "Tstat",
                 "var", "p.value", "Tstat_c", "p.value_c", "varT_c",
                 "BETA_c", "SE_c") 
  
  # subset columns
  d <- d[,colnames(d) %in% keep_cols, with = FALSE]
  marker <- apply( d[ ,c("CHR", "POS", "Allele1", "Allele2") ] , 1 , paste , collapse = ":" )
  rsid <- d$MarkerID
  pvalue <- d[[col_pvalue]]
  tstat <- d[[col_tstat]]
  zstat <- d[[col_beta]] / d[[col_se]]
  logp <- z_to_log10p(zstat) 
  p <- 2 * pnorm(abs(zstat), lower.tail = FALSE)

  # setup new data.table with columns
  dnew <- data.table(
    cur_marker=marker,  
    cur_rsid=rsid, 
    cur_pvalue=pvalue, 
    cur_tstat=tstat, 
    cur_zstat=zstat,
    cur_log10p=logp,
    cur_zstat_p=p
  )

  # keep original
  d <- cbind(dnew, d)

  # perform rows based on markers that pass thresholds
  p_cutoff <- as.numeric(args$p_cutoff)
  bool <- d$cur_pvalue < p_cutoff
  d <- d[bool, ]
  
  # check for any INF estimates
  bool_inf <- is.infinite(d$cur_zstat)
  n_inf <- sum(bool_inf)
  d <- d[!bool_inf, ]
  msg <- paste(n_inf, "rows are have infinite beta/SEs, these have been discarded.")
  if (n_inf > 0) write(msg, stderr())

  # sort by most significant at bottom. Note that we are sorting by t-stats
  # since floating point precision is not good enough for SAIGE to do extreme P-values
  #d <- d[rev(order(abs(d$pvalue))), ]
  #d <- d[order(abs(d$cur_tstat)), ]
  d <- d[order(abs(d$cur_zstat)), ]
  #d <- d[order(abs(d$cur_log10p)), ]

  d_subset <- tail(d, 1)
  msg <- paste("Lowest P-value/log10(p-value)/tstat detected was:", d_subset$cur_pvalue,"/", d_subset$cur_log10p ,"/", d_subset$cur_tstat)
  write(msg, stderr())


  write(paste("Writing", args$out_file), stderr())
  fwrite(d, args$out_file, sep = "\t", col.names = TRUE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--spa_file", default=NULL, help = "?")
parser$add_argument("--out_file", default=NULL, help = "?")
parser$add_argument("--p_cutoff", default=NULL, help = "?")
parser$add_argument("--col_pvalue", default="p.value", help = "?")
parser$add_argument("--col_pvalue_cond", default="p.value_c", help = "?")
args <- parser$parse_args()

main(args)

