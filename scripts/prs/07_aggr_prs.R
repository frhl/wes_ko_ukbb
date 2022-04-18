
devtools::load_all('utils/modules/R/prstools')
library(argparse)

main <- function(args){

  print(args)
  stopifnot(dir.exists(args$in_dir))
  stopifnot(file.exists(args$path_phenos)) 
  
  phenos <- unlist(strsplit(readLines(args$path_phenos), split = "\\s+", perl = TRUE)) 
  files <- list.files(args$in_dir, pattern = '.txt.gz', full.names = TRUE)
  print(phenos)
  print(files)

   pheno_list <- lapply(phenos, function(pheno){
      
      if (is.null(args$grep)){
         regex <- paste0('_',pheno,'_')
         files_regex <- files[grepl(regex, files)]
      } else {
         files_regex <- files[grepl(args$grep, files)]
      }
         
      if (length(files_regex) > 0) {
          
          write(paste0('Aggregating ',pheno,'..'),stdout())
          # append into a single 2 x n matrix
          d <- do.call(rbind, lapply(1:22, function(x) {
              f <- files_regex[grepl(paste0("chr",x,".txt.gz"),files_regex)]
              if (length(f)){
                  d <- fread(f)
                  d$chr <- paste0("chr",x)
                  return(d)
              } else {
                write(paste0("Warning: ", pheno, "chr", x," does not exists!"), stderr())
              }
          }))

          # Aggregate PRS into a n x 22 matrix (n = samples)
          M <- dcast(sid ~ chr, value.var = "prs", data = d)

          # sum up PRS score
          M[is.na(M)] <- 0
          M <- data.table(
              sid = M$sid, 
              pgs = rowSums(M[,2:ncol(M)]),
              phenotype = pheno
          )
          
          return(M)
          
      } else {
          
          return(NULL)
      }
      
  })

  # combine into single matrix
  d <- do.call(rbind, pheno_list)
  M <- dcast(sid ~ phenotype, value.var = "pgs", data = d)
  fwrite(M, paste0(args$out_prefix, '.txt.gz'), sep = '\t') 

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Directory for (chromosome-wise) PRS")
parser$add_argument("--path_phenos", default=NULL, required = TRUE, help = "Path to phenotype file with phenotype names as header")
parser$add_argument("--grep", default=NULL, help = "Subset in_dir files by command")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









