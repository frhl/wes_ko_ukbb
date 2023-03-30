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
  stopifnot(dir.exists(args$ld_dir))
  stopifnot(args$chrom %in% paste0("chr",1:22))
  stopifnot(dir.exists(dirname(args$out_prefix)))
  stopifnot(dir.exists(dirname(args$path_betas)))

  tic(paste0("LDPred2 ", basename(args$out_prefix)))
  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)

  # laod ldsc and qced GWAS
  ldsc <- readRDS(args$ldsc)
  h2 <- ldsc$coefficients$estimate[2]
  pvalue <- ldsc$coefficients$pvalue[2]
  n_eff <- min(unique(ldsc$gwas$n_eff))
  gwas <- ldsc$gwas

  # check estimates and ensure that
  # the trait is actually heritable 
  #stopifnot(pvalue < 1e-5)
  stopifnot(!is.null(gwas)) 

  # load prediction file (note, that this needs to be done before loading
  # LD-matrix, so that it can be subsetted accordinly).
  pred <- load_bigsnp_from_bed(args$pred)
  gwas <- gwas[gwas$rsid %in% pred$map$rsid,]

  # need to keep track of matching SNPs.
  n_pred_in_gwas <- sum(pred$map$rsid %in% gwas$rsid)
  total <- nrow(pred$map)
  pct <- paste0("(",round((n_pred_in_gwas/total)*100, 2),"%)")
  msg <- paste("Only",n_pred_in_gwas,"of",total,"SNPs", pct,"GWAS and prediction SNVs:", args$out_prefix)
  if (n_pred_in_gwas < 1000) stop(msg)
      
  # get SNP correlations and LD
  snp <- get_single_ld_matrix(gwas, chr = args$chrom, ld_dir = args$ld_dir)
 
  # match GWAS with snp-correlation map
  gwas <- gwas[na.omit(match(snp$map$marker, gwas$marker)), ]

  # check that LD-matrix markers and gwas markers have overlap
  # Check that ordering of markers are actually matching
  stopifnot(all(gwas$marker %in% snp$map$marker))
  stopifnot(all(snp$map$marker %in% gwas$marker)) 
  stopifnot(sum(gwas$marker == snp$map$marker) / nrow(gwas) == 1)
  
  # load data to be used for prediction 
  bfile <- tempfile(tmpdir = dirname(args$out_prefix))
  pred <- match_bigsnp_with_gwas(obj=pred, gwas=gwas, bfile=args$tmp_bfile)
  genotypes <- pred$genotypes  
  indicies <- pred$gwas_indicies

  # check that we have genotypes
  stopifnot(!is.null(genotypes))
  stopifnot(!is.null(indicies))
    
  # check that things are mathced
  stopifnot(pred$map$rsid == gwas$rsid)

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

  # save chains
  multi_auto <- readRDS(path_betas)
  final_beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
  final_beta_auto <- gwas$final_beta_auto

  # perform prediction
  final_pred <- big_prodVec(
     genotypes, 
     final_beta_auto,
     center = means,
     scale = sds,
     ncores = NCORES)

  # count PRS that could not be generated
  missing_prs <- round(sum(is.na(final_pred))/length(final_pred)*100, 2)

  # check if all PRS values are predicted as zero
  if (all(is.na(final_pred))) stop(paste("Error! Predicted PRS scores are all NA for", args$out_prefix ))
  if (all(as.numeric(na.omit(final_pred))==0)) stop(paste("Error! Predicted PRS scores are all zero for", args$out_prefix ))

  # save results
  PGS <- data.table(
    sid = pred$sid, 
    prs = final_pred
  )

  write(paste0(args$pred, ".. done! Writing to ", args$out_prefix, ".txt.gz"), stdout())
  fwrite(PGS, file = paste0(args$out_prefix,".txt.gz"), sep = '\t')
  fwrite(model, file = paste0(args$out_prefix,".model"), sep = '\t')
  toc()
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "chromosome input")
parser$add_argument("--pred", default=NULL, required = TRUE, help = "Path to plink (bed) for PGS prediction")
parser$add_argument("--standardized_gt", default=1, required = FALSE, help = "Should genotypes be standardized?")
parser$add_argument("--tmp_bfile", default=NULL, required = TRUE, help = "File path to temporary backing files")
parser$add_argument("--ld_dir", default=NULL, required = TRUE, help = "Path to directory with pre-calcualted SNP correlations and LD (.rds files)")
parser$add_argument("--impute", default=NULL, required = TRUE, help = "Should missing genotypes be imputed? (See https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--path_betas", default=NULL, required = TRUE, help = "What is the path to the betas that should be used?")
args <- parser$parse_args()

main(args)



