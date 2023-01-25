#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

    write("Reading phenotypes..", stderr())
    stopifnot(file.exists(args$phenos))
    d <- fread(args$phenos)
    
    d$eid <- as.character(d$eid)
    header <- readLines(args$header)
    header <- header[!grepl("primary", header)]
    stopifnot(length(header) < 1000)

    # get unique combinations
    header_mat <- data.table::CJ(header, header, unique = TRUE)
    rownames(header_mat) <- NULL
    colnames(header_mat) <- c("Var1", "Var2")

    # setup comparison
    tested <- c()
    indexes <- 1:nrow(header_mat)
    tab <- do.call(rbind, lapply(indexes, function(idx){
        x <- as.character(header_mat$Var1[idx])
        y <- as.character(header_mat$Var2[idx])
        target <- paste0(sort(c(x,y)), collapse = "-")
        if (!(target %in% tested)) {
            write(target, stderr())
            tested <- c(tested, target)
            stopifnot(x %in% header)
            stopifnot(y %in% header)
            dx <- data.table(eid=d$eid, p=d[[x]])
            dy <- data.table(eid=d$eid, p=d[[y]])
            dx <- dx[!is.na(dx$p),]
            dy <- dy[!is.na(dy$p),]
            
            dx_cases <- dx$eid[(as.logical(dx$p) == TRUE)]
            dx_ctrls <- dx$eid[(as.logical(dx$p) == FALSE)]
            n_dx_cases <- length(dx_cases)
            n_dx_ctrls <- length(dx_ctrls)

            dy_cases <- dy$eid[(as.logical(dy$p) == TRUE)]
            dy_ctrls <- dy$eid[(as.logical(dy$p) == FALSE)]
            n_dy_cases <- length(dy_cases)
            n_dy_ctrls <- length(dy_ctrls)

            overlap_cases <- intersect(dx_cases, dy_cases)
            overlap_ctrls <- intersect(dx_ctrls, dy_ctrls)        
            n_overlap_cases <- length(overlap_cases)
            n_overlap_ctrls <- length(overlap_ctrls)

            x_diff_y <- setdiff(dx_cases, dy_cases)
            y_diff_x <- setdiff(dy_cases, dx_cases)
            n_x_diff_y <- length(x_diff_y)
            n_y_diff_x <- length(y_diff_x)
            
            out <- data.table(
              x = x,
              y = y,
              defined_x = n_dx_cases + n_dx_ctrls,
              defined_y = n_dy_cases + n_dy_ctrls,
              cases_x = n_dx_cases,
              cases_y = n_dy_cases,
              cases_overlap_xy = n_overlap_cases,
              ctrls_overlap_xy = n_overlap_ctrls,
              cases_x_diff_y = n_x_diff_y,
              cases_y_diff_x = n_x_diff_y,
              cases_x_diff_y_div_x = n_x_diff_y/n_dx_cases,
              cases_y_diff_x_div_y = n_y_diff_x/n_dy_cases
            )
            return(out)
        }
    }))

    outfile <- paste0(args$out_prefix,".txt.gz")
    write(paste0("writing ", outfile,".."), stderr())
    fwrite(tab, outfile)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenos", default=NULL, help = "chromosome")
parser$add_argument("--header", default=NULL, help = "chromosome")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

