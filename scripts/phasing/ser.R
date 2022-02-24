devtools::load_all("utils/modules/R/phasingtools")

main <- function(args){

  print(args)
  stopifnot(dir.exists(args$in_dir))
  stopifnot(dir.exists(args$out_prefix))
  
  files <- list.files(args$in_dir, pattern = args$in_ext, full.names = TRUE)
  if (!is.null(grep)) files <- files[grepl(args$grep, basename(files))]
  stopifnot(length(files) > 0)

  # simple analysis of non-hacked BCFtools
  if (!is.null(args$simple)){
  
    M <- do.call(rbind, lapply(files, function(f){
        d <- fread(f) 
        d <- summarize_bcftools_trio_stats(d)  
        d$chr <- str_extract(f, 'chr[0-9]+')
        return(d)
    }))

    return(M)
  }

  # Site-specific analysis
  


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--in_ext", default='*.txt', help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--grep", default=NULL, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--simple", default=NULL, action = 'store_true', help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


