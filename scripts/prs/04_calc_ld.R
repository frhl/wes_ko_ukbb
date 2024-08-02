
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html
# issues with parallel: https://github.com/privefl/bigstatsr/issues/90

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)


main <- function(args){

  print(args)
  stopifnot(file.exists(args$bed))

  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)
  
  # load data required for setting up LD-matrix 
  ld_data <- load_bigsnp_from_bed(args$bed)
  POS2 = ld_data$POS2
  G = ld_data$G

  chrs <- paste0("chr",1:22)
  for (chr in chrs) {
      write(paste('Calculating LD for',chr,'..'),stderr())
      ind.chr <- which(ld_data$map$chr == chr)
      stopifnot(length(ind.chr) > 0)
      corr0 <- snp_cor(
              G,
              ind.col = ind.chr,
              ncores = NCORES,
              infos.pos = POS2[ind.chr],
              size = 3/1000
          )
      if (!is.null(args$compress)){
        saveRDS(corr0, paste0(args$out_prefix, "_",chr,".rda"), compress = args$compress)
      } else {
        saveRDS(corr0, paste0(args$out_prefix, "_",chr,".rds"))
      }
      
      fwrite(ld_data$map[ind.chr, ], paste0(args$out_prefix, "_",chr,'.txt.gz'))

      #if (chr == "chr1") {
      #    ld <- Matrix::colSums(corr0^2)
      #    corr <- as_SFBM(corr0, args$out_prefix)
      #} else {
      #    ld <- c(ld, Matrix::colSums(corr0^2))
      #    corr$add_columns(corr0, nrow(corr))
      #}
  }  
 
  snp_corr <- (list(ld = ld, corr = corr))

  # save with link to sfbm file
  #saveRDS(snp_corr, file = paste0(args$out_prefix, ".rda"), compress = "xz")
  #fwrite(ld_data$map, paste0(args$out_prefix, ".txt.gz"))

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--bed", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--compress", default=NULL, help = "what type of compression should be applied?")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


