# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)
library(ggplot2)

main <- function(args){

    stopifnot(file.exists(args$ldsc))
    stopifnot(file.exists(args$path_betas))
    stopifnot(file.exists(args$path_betas_map))
    stopifnot(dir.exists(dirname(args$out_prefix)))

    # laod ldsc and qced GWAS
    ldsc <- readRDS(args$ldsc)
    gwas <- ldsc$gwas

    # load betas and do chain QC
    multi_auto <- readRDS(args$path_betas)
    (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
    (keep <- (range > (0.95 * quantile(range, 0.95, na.rm=TRUE))))
    keep[is.na(keep)] <- FALSE
    stopifnot(sum(keep) > 0)

    # obtain weights
    betas <- fread(args$path_betas_map)
    final_beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    betas$beta <- final_beta_auto

    # write combined QC-ed PRS weights
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(betas, outfile, sep="\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_betas", default=NULL, required = TRUE, help = "Path to LDPred2 chains")
parser$add_argument("--path_betas_map", default=NULL, required = TRUE, help = "Path to SNPs corresponding to weight indicies")
parser$add_argument("--ldsc", default=NULL, required = TRUE, help = ".rds object containing QCed GWAS and ldsc heritability estimates")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)



