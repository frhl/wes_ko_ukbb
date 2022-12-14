
library(argparse)
library(data.table)
library(bigstatsr)

main <- function(args){

    print(args)
    stopifnot(file.exists(args$phenotype_bin)) 
    stopifnot(dir.exists(args$directory))
    stopifnot(dir.exists(dirname(args$out_prefix)))
 
    # Load polygenic risk scores
    files <- list.files(args$directory, pattern = ".txt.gz", full.names = TRUE)
    files <- files[!grepl("chrom",files)]

    # combine files
    lst <- lapply(files, fread)
    mrg <- Reduce(merge, lst)
    mrg$eid <- mrg$sid

    # read binary phenotypes
    d_bin <- fread(args$phenotype_bin)
    d_bin <- d_bin[d_bin$eid %in% mrg$eid,]
    cols_bin <- colnames(d_bin)

    # subset cols to binary
    cols <- gsub("_pgs","",colnames(mrg))
    cols <- cols[cols %in% cols_bin]
    d_bin <- d_bin[,colnames(d_bin) %in% cols, with = FALSE]
    dt <- merge(mrg, d_bin, by = 'eid')

    # calculate AUC 
    cols <- cols[!cols %in% "eid"]
    lst <- lapply(cols, function(col){
        write(paste("calculating auc for", col), stderr())
        col_pgs <- paste0(col,'_pgs')
        boot <- dt[,colnames(dt) %in% c('eid',col,col_pgs), with = FALSE]
        boot <- boot[!is.na(boot[[col]]) & !is.na(boot[[col_pgs]]),]
        auc <- AUCBoot(
            pred = boot[[col_pgs]],
            target = as.numeric(boot[[col]]),
            nboot = 10000,
            seed = 1995
        )
        
        cases <- sum(boot[[col]]==1)
        controls <- sum(boot[[col]]==0)
        boot <- data.table(t(auc))
        colnames(boot) <- paste0("auc_",tolower(colnames(boot)))
        colnames(boot) <- gsub("\\%","_pct", colnames(boot))
        colnames(boot) <- gsub("\\.","_", colnames(boot))
        boot$pred_cases <- cases
        boot$pred_controls <- controls
        boot$pred_n <- cases + controls
        boot$phenotype <- col
        
        return(boot)
        
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
parser$add_argument("--phenotype_bin", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









