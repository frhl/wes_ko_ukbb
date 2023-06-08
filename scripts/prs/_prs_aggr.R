
devtools::load_all('utils/modules/R/prstools')
library(argparse)

main <- function(args){

  print(args)
  stopifnot(dir.exists(args$in_dir))
  stopifnot(dir.exists(args$out_dir))
  
  print(paste("Looking for chromosomes in", args$in_dir))
  files <- list.files(args$in_dir, pattern = paste0("pgs.",args$phenotype,".chr[0-9]+.txt.gz"), full.names = TRUE)
  
  n <- length(files)
  if (n < 22) {
    stop(paste("Found ",n, "files for", args$phenotype))
  }

  # combine into pgs matrix (chr
  d <- do.call(rbind, lapply(1:22, function(x) {
      f <- files[grepl(paste0("\\.chr",x,"\\.txt\\.gz"),files)]
      write(paste("Attempting to open", f), stderr())
      if (length(f)){
          d <- fread(f)
          d$chr <- paste0("chr",x)
          return(d)
      } else {
        write(paste0("Warning: ", args$phenotype, "chr", x," does not exists!"), stderr())
      }
  }))

  # Aggregate PRS into a n x 22 matrix (n = samples)
  M <- dcast(sid ~ chr, value.var = "prs", data = d)
  M <- M[,c('sid', paste0("chr",1:22))]

  # sum up PRS score
  M[is.na(M)] <- 0
 
  # write matrix with each chromosome 
  outfile1 <- paste0(args$out_dir, "/", args$phenotype, "_pgs_chrom.txt.gz")
  write(paste("writing to", outfile1), stdout())
  fwrite(M, outfile1, sep = "\t") 
 
  # write combined matrix
  M <- data.table(
    sid = M$sid,
    pgs = rowSums(M[,2:ncol(M)]))
  colnames(M)[-1] <- paste0(args$phenotype, '_', colnames(M)[-1])
  outfile2 <- paste0(args$out_dir, "/", args$phenotype, "_pgs.txt.gz")
  fwrite(M, outfile2, sep = "\t")



}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Directory for (chromosome-wise) PRS")
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "")
parser$add_argument("--regex", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_dir", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









