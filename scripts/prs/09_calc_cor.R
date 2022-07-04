
library(argparse)
library(data.table)
library(bigstatsr)

main <- function(args){

    print(args)
    stopifnot(file.exists(args$phenotype_cts)) 
    stopifnot(dir.exists(args$directory))
    stopifnot(dir.exists(dirname(args$out_prefix)))
 
    # Load polygenic risk scores
    files <- list.files(args$directory, pattern = ".txt.gz", full.names = TRUE)
    files <- files[!grepl("chrom",files)]

    # comctse files
    lst <- lapply(files, fread)
    mrg <- Reduce(merge, lst)
    mrg$eid <- mrg$sid

    # read ctsary phenotypes
    d_cts <- fread(args$phenotype_cts)
    d_cts <- d_cts[d_cts$eid %in% mrg$eid,]
    cols_cts <- colnames(d_cts)

    # subset cols to ctsary
    cols <- gsub("_pgs","",colnames(mrg))
    cols <- cols[cols %in% cols_cts]
    d_cts <- d_cts[,colnames(d_cts) %in% cols, with = FALSE]
    dt <- merge(mrg, d_cts, by = 'eid')

    # calculate AUC 
    cols <- cols[!cols %in% "eid"]
    lst <- lapply(cols, function(col){
        print(col)
        col_pgs <- paste0(col,'_pgs')
        cur_dt <- dt[,colnames(dt) %in% c('eid',col,col_pgs), with = FALSE]
        cur_dt <- cur_dt[!is.na(cur_dt[[col]]) & !is.na(cur_dt[[col_pgs]]),]
        f <- as.formula(paste0(col, "~", col_pgs))
        fit <- summary(lm(f, data = cur_dt))
        correlation <- cor(cur_dt[[col]], cur_dt[[col_pgs]])
        
        # perform Z-test
        estimate <- fit$coefficients[2,1]
        stderr <- fit$coefficients[2,2]
        zscore <- estimate / stderr
        pvalue <- 2 * pnorm(abs(zscore), lower.tail = FALSE, log.p = FALSE)
        log_pvalue <- 2 * pnorm(abs(zscore), lower.tail = FALSE, log.p = TRUE)
        d_out <- data.frame(
            phenotype = col,
            correlation = correlation,
            estimate = estimate,
            stderr = stderr,
            p = pvalue,
            log_p = log_pvalue,
            pred_n = nrow(d_cur)
        )
        return(d_out) 
    })

    # write final file
    final <- data.table(do.call(rbind, lst))
    outfile <- paste0(args$out_prefix, ".txt.gz")
    write(paste("writing",outfile),stdout())
    fwrite(final, outfile, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--directory", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--phenotype_cts", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









