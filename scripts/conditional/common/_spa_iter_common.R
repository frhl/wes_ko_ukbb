#!/usr/bin/env Rscript

library(argparse)
library(data.table)

z_to_log10p <- function(z) {
      abs((pnorm(-abs(z),log.p=TRUE)+log(2))/log(10))
}

main <- function(args){

    # read and check integratity of file
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
    d_marker <- d[ ,c("CHR", "POS", "Allele1", "Allele2") ]

    # create new data.table containing currently used columns
    marker <- apply(d_marker, 1, paste, collapse = ":" )
    marker <- gsub("\\ +", "", marker)
    dt <- data.table(cur_marker = marker)
    dt$cur_pvalue <- d[[col_pvalue]]
    dt$cur_tstat <- d[[col_tstat]]
    dt$cur_zstat <- d[[col_beta]] / d[[col_se]]
    dt$cur_zstat_pvalue <- 2 * pnorm(abs(dt$cur_zstat), lower.tail = FALSE)
    dt$cur_zstat_log10_pvalue <- z_to_log10p(dt$cur_zstat)

    # do some checks
    #stopifnot(round(cor(dt$cur_pvalue, dt$cur_zstat_pvalue), 6) == 1)
    marker_zstat_p <- dt$cur_marker[which(dt$cur_zstat_log10_pvalue == max(dt$cur_zstat_log10_pvalue))]
    marker_cur_p <- dt$cur_marker[which(dt$cur_pvalue == min(dt$cur_pvalue))]

    # lowest P
    min_p <- min(dt$cur_pvalue)

    # combine with original
    d <- data.table(cbind(dt, d))
    p_cutoff <- as.numeric(args$p_cutoff)
    bool_keep <- d$cur_pvalue < p_cutoff
    d <- d[bool_keep, ]
    finished <- nrow(d) == 0
    if (finished) write(paste0("SUCCESS! No marker markers left to condition on (Current min P-value = ",min_p,")"), stderr())


    # ensure that most significant marker is at bottom
    #d <- d[rev(order(d$cur_zstat_pvalue)),]
    write("Sorting file by Tstat (Most sig. at bottom)", stderr())
    d <- d[order(abs(d$cur_tstat)),]


    use_col_names = nrow(d) > 0
    write(paste("Writing", args$out_file), stderr())
    fwrite(d, args$out_file, sep = "\t", col.names = use_col_names)

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

