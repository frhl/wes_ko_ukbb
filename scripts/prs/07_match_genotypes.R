# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)
library(ggplot2)
library(tictoc)

main <- function(args){

  stopifnot(file.exists(args$pred))
  stopifnot(!file.exists(args$tmp_bfile))
  stopifnot(args$chrom %in% paste0("chr",1:22))
  stopifnot(dir.exists(dirname(args$out_prefix)))
  stopifnot(file.exists(args$path_weights))
  
  # setup parallel environment
  tic(paste0("LDPred2 ", basename(args$out_prefix)))
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)

  # get weights
  weights <- fread(args$path_weights)
  pred <- load_bigsnp_from_bed(args$pred)
  weights <- weights[weights$marker %in% pred$map$rsid, ]

  # check that pred file hsa enough markers
  n_pred_in_weights <- sum(pred$map$rsid %in% weights$marker)
  total <- nrow(pred$map)
  pct <- paste0("(",round((n_pred_in_weights/total)*100, 2),"%)")
  msg <- paste("Only",n_pred_in_weights,"of",total,"SNPs", pct,"GWAS and prediction SNVs:", args$out_prefix)
  if (n_pred_in_weights < 1000) stop(msg)
  
  # load data to be used for prediction 
  bfile <- tempfile(tmpdir = dirname(args$out_prefix))
  pred <- match_bigsnp_with_gwas(obj=pred, gwas=weights, bfile=args$out_prefix)
  genotypes <- pred$genotypes  
  indicies <- pred$gwas_indicies

  # check that we have genotypes
  stopifnot(!is.null(genotypes))
  stopifnot(!is.null(indicies))
    
  # check that things are mathced
  stopifnot(pred$map$rsid == weights$marker)

  # dimensions  
  cols <- genotypes$`.->ncol` # variants
  rows <- genotypes$`.->nrow` # samples
  missing_gt <- NA

  print(paste(cols, "variants"))

  # check if any SNVs are invariant
  snv_var <- as.numeric(unlist(lapply(1:cols, function(i) var(genotypes[,i], na.rm = TRUE))))
  invariant_index <- which(snv_var == 0)                                    
  if (length(invariant_index) > 0) stop(paste(length(invariant_index),"SNVs are invariant!"))

  # need to impute missing SNPs
  if (!is.null(args$impute)){
    write(paste0("Imputing genotypes in ", args$out_prefix), stderr())
    sum_cols <- lapply(1:cols, function(i) return(sum(is.na(genotypes[,i])))) 
    missing_gt <- sum(unlist(sum_cols))
    genotypes <- snp_fastImputeSimple(genotypes, method = args$impute) 
  }

  # standardize genotypes 
  means <- NULL; sds <- NULL
  if (args$standardized_gt){
     write(paste0("Standardizing genotypes in ", args$out_prefix), stderr())
     means <- as.numeric(unlist(lapply(1:cols, function(i) mean(genotypes[,i], na.rm = TRUE))))
     sds <- as.numeric(unlist(lapply(1:cols, function(i) sd(genotypes[,i], na.rm = TRUE))))
     stopifnot(sum(is.na(means))==0)
     stopifnot(sum(is.na(sds))==0)
  }
 

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "chromosome input")
parser$add_argument("--pred", default=NULL, required = TRUE, help = "Path to plink (bed) for PGS prediction")
parser$add_argument("--standardized_gt", default=1, required = FALSE, help = "Should genotypes be standardized?")
parser$add_argument("--impute", default=NULL, required = TRUE, help = "Should missing genotypes be imputed? (See https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--path_weights", default=NULL, required = TRUE, help = "What is the path to the betas that should be used?")
args <- parser$parse_args()

main(args)



