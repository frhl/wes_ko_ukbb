
library(argparse)
library(data.table)
library(ggplot2)
library(boot)

# a function to calculate non-parametric botostrap correlation
boot_cor_ci <- function(data, col1, col2, n = 1000, method='pearson'){
    data <- data.frame(data)
    b <- boot(data, 
      statistic = function(data, i) {
        cor(data[i, col1], data[i, col2], method='pearson')
      },
      R = n
    )
    out <- data.frame(t(matrix(quantile(b$t, probs = c(0.025, 0.975)))))
    colnames(out) <- c("ci0025","ci975")
    out$mean <- mean(b$t)
    out$col1 <- col1
    out$col2 <- col2
    return(out)
}

# a function to perform a standard z-test
ztest <- function(d, a1, a2){
    f <- as.formula(paste0(a1, "~", a2))
    fit <- lm(f, data = d)
    coef <- coefficients(summary(fit))
    est <- coef[2,1]
    stderr <- coef[2,2]
    z <- est / stderr
    pvalue <- 2 * pnorm(abs(z), lower.tail = FALSE)
    return(pvalue)
}


main <- function(args){

    print(args)
    stopifnot(file.exists(args$pgs_dir)) 
    stopifnot(dir.exists(dirname(args$out_dir)))
    autosomes <- paste0("chr",1:22)
    files <- list.files(args$pgs_dir, pattern = "chrom\\.txt\\.gz", full.names = TRUE)
    outpdf <- paste0(args$out_dir,"_qc.pdf")
    pdf(outpdf, width = 10, height = 8)
    for (f in files){
        write(paste("running", f), stderr())        
        d <- fread(f) 
        # setup files
        tested <- c()
        M <-do.call(rbind, lapply(autosomes, function(a1){
            do.call(rbind, lapply(autosomes, function(a2){
                    testing <- paste0(sort(c(a1, a2)))
                    out <- NULL
                    if ((a1 != a2) & (!testing %in% tested)) {out <- data.frame(boot_cor_ci(d, a1, a2), pvalue = ztest(d, a1, a2))}
                    tested <- c(tested, testing)
                    return(out)
            }))
        }))

        # prepare writing file
        M <- data.table(M)
        M$filepath <- basename(f)
        f_wo_ext <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(f))) 
        outfile <- file.path(args$out_dir, paste0(f_wo_ext,".txt.gz"))
        fwrite(M, outfile, sep = "\t")
    }
    graphics.off()

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--pgs_dir", default=NULL, required = TRUE, help = "polygenic risk score directory")
parser$add_argument("--out_dir", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









