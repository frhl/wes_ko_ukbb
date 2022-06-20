# investigate the correlation of PRS between samples


devtools::load_all('utils/modules/R/prstools')
library(argparse)
library(ggplot2)

main <- function(args){

  print(args)
  stopifnot(file.exists(args$prs)) 
  stopifnot(file.exists(args$phenotypes)) 

  pgs <- fread(args$prs)
  if (args$sample_unif > 0){
    set.seed(args$seed)
    n <- nrow(pgs)
    stopifnot(args$sample_unf < n)
    sample_idx <- sample(1:n, args$sample_unif, replace  = FALSE)
    pgs <- pgs[sample_idx, ]
  } 

  
  form <- paste("~",paste0(paste0("chr", 1:22), collapse = "+"), collapse = "")
  p <- pairs(as.formula(form), data = pgs)
  
  
  #plt <- ggplot(d, aes(x=pgs, fill=is_case)) +
  #  geom_histogram(color = 'black') +
  #  labs(fill="Case") + 
   # xlab('Polygenic Risk Score') +
   # ylab("Frequency") +
   # theme_bw() +
   # facet_wrap(~phenotype)  

  ggsave(paste0(args$prefix, ".pdf"), plt, width = 12, height = 14)
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--prs", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--phenotypes", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--sample_unif", default=0, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--seed", default=42, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









