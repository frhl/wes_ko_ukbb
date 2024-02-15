#!/usr/bin/env Rscript

library(argparse)
library(data.table)

r_temp = "data/tmp/rtmp"

main <- function(args){
    
    stopifnot(file.exists(args$input_path))
    stopifnot(dir.exists(dirname(args$out_prefix)))
    #stopifnot(is.numeric(as.numeric(args$true_tstat)))
    stopifnot(is.numeric(as.numeric(args$true_p)))
    
    # get true values from non-permuted data
    #true_t <- as.numeric(args$true_tstat)
    #true_p <- as.numeric(args$true_p)
    #true_AC <- as.numeric(args$true_AC)
    marker_id <- as.character(args$marker_id)

    # read input
    d <- fread(args$input_path)

    # We include the "real" marker in the VCF to extract 
    # the reference/true P-value and T-statistic
    true_p_marker_name <- "actual"
    stopifnot(true_p_marker_name %in% d$MarkerID)
    d_actual <- d[d$MarkerID == true_p_marker_name,]
    if ("p.value_c" %in% colnames(d_actual)){ 
      true_p <- as.numeric(d_actual$p.value_c)
      true_t <- as.numeric(d_actual$Tstat_c)
      true_v <- as.numeric(d_actual$var_c)
    } else {
      true_p <- as.numeric(d_actual$p.value)
      true_t <- as.numeric(d_actual$Tstat)
      true_v <- as.numeric(d_actual$var)
    }

    # if a merged VCF with two actual/real
    # markers for the non-permuted P is read 
    # we need to take the unique
    true_p <- unique(true_p)
    true_t <- unique(true_t)
    true_v <- unique(true_v)
    stopifnot(length(true_p) == 1)
    stopifnot(length(true_t) == 1)
    stopifnot(length(true_v) == 1)

    # exclude real markers (non-permuted stuff)
    stopifnot("MarkerID" %in% colnames(d))
    bool_real_marker <- !grepl("ENSG", d$MarkerID)
    n_real_marker <- sum(bool_real_marker)
    if (n_real_marker > 0) {
        write(paste("Excluded", n_real_marker, "real marker(s)."), stderr()) 
        d <- d[!bool_real_marker, ]
    }


    # get the p-values for the resulting 
    # permuted markers
    out_marker <- unique(d$MarkerID)
    stopifnot(length(out_marker) == 1)
    if ("p.value_c" %in% colnames(d)){ 
      pvalue <- as.numeric(d$p.value_c)
      tstat <- as.numeric(d$Tstat_c)
      var <- as.numeric(d$var_c)
    } else {
      pvalue <- as.numeric(d$p.value)
      tstat <- as.numeric(d$Tstat)
      var <- as.numeric(d$var)
    }

    # data.frame of true P-values and permuted P-values
    dt_true <- data.table(Tstat=true_t, p=true_p, var=true_v,  out_marker=out_marker, is_permuted=0)
    dt_perm <- data.table(Tstat=tstat, p=pvalue,  var=var, out_marker=out_marker, is_permuted=1)
    dt <- rbind(dt_true, dt_perm)
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(dt, outfile, sep = "\t")
    
    # count how many times the original P-value could be found
    idx_same_p <- which(true_p %in% dt$p[dt$is_permuted == 1])  
    n_same_p <- length(idx_same_p)
    idx_string <- paste0(idx_same_p, collapse = ",")

    if (args$alternative == "two.sided"){
        count_p_ge_true <- sum(true_p >= pvalue)/length(pvalue) 
    } else if (args$alternative == "greater"){
        count_p_ge_true <- sum(true_t >= tstat)/length(tstat) 
    } else {
        stop("arg 'alternative' must be either 'two.sided' or 'greater'")
    }
    
    empirical_p <- count_p_ge_true
    write(empirical_p, stdout())


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "path to saige (merge) file of markers")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for output")
parser$add_argument("--alternative", default="two.sided", help = "Character string specifying the alternative hypothesis, must be one of 'two.sided' (default) or 'greater'.")
parser$add_argument("--marker_id", default=NULL, help = "the true marker name")
parser$add_argument("--exclude_real_markers", default=FALSE, action="store_true", help = "Exclude real conditioning markers (from conditioning analysis)")
args <- parser$parse_args()

main(args)

