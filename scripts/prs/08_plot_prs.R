
devtools::load_all('utils/modules/R/prstools')
library(argparse)
library(ggplot2)

main <- function(args){

  print(args)
  stopifnot(file.exists(args$prs)) 
  stopifnot(file.exists(args$phenotypes)) 

  pgs <- fread(args$prs)
  pgs_phenotypes <- colnames(pgs)[-1]
  phenotypes <- fread(args$phenotypes)

  d <- do.call(rbind,lapply(pgs_phenotypes, function(phenotype){
    df <- phenotypes[,c('eid',phenotype), with = FALSE]
    colnames(df) <- c('sid','phenotype')
    mrg <- merge(pgs, df, all.x = TRUE)
    colnames(mrg) <- c('sid','pgs','is_case')
    mrg$phenotype <- phenotype
    return(mrg)
  }))
  
  plt <- ggplot(d, aes(x=pgs, fill=is_case)) +
    geom_histogram(color = 'black') +
    labs(fill="Case") + 
    xlab('Polygenic Risk Score') +
    ylab("Frequency") +
    theme_bw() +
    facet_wrap(~phenotype)  

  ggsave(paste0(args$prefix, ".pdf"), plt, width = 12, height = 14)
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--prs", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--phenotypes", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









