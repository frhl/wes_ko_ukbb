
devtools::load_all('utils/modules/R/prstools')
library(argparse)

main <- function(args){

  print(args)

  # get the right files
  files <- list.files(args$in_dir, pattern = '.txt.gz', full.names = TRUE)
  files <- files[grepl(files, args$grep)]
  
  # append into a single 2 x n matrix
  d <- do.call(rbind, lapply(1:22, function(x) {
      f <- files[grepl(paste0("chr",x,".txt.gz"),files)]
      print(f)
      if (length(f)){
          d <- fread(f)
          d$chr <- paste0("chr",x)
          return(d)
      }
  }))
  
  # Aggregate PRS into a n x 22 matrix (n = samples)
  M <- dcast(sid ~ chr, value.var = "prs", data = d)
  
  # sum up PRS score 
  M[is.na(M)] <- 0
  M <- data.table(sid = M$sid, pgs = rowSums(M[,2:ncol(M)]))
  
  # write file
  fwrite(M, paste0(args$out_prefix, '.txt.gz'), sep = '\t') 

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--grep", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









