#' @title read hail summary statistics
#' @param path string file to summary statistics generated with HAIL
#' @param remove_failed boolean. remove SNPs with NAs in betas?
#' @return data.table of ldpred2 ready summary statistics
#' @export

read_hail_sumstat <- function(path, trait, remove_failed = TRUE){
  # Read summary statistics and format to reference
  # from hail to ldpred2 matching coluns.
  stopifnot(file.exists(path))
  stopifnot(trait %in% c("binary", "cts"))

  sumstats <- bigreadr::fread2(path)
  colnames(sumstats)[colnames(sumstats)=="n_total"] <- 'n' # tmp fix
  # make artificial rsid
  sumstats$rsid <-
      apply(
          data.frame(
              locus = sumstats$locus,
              a0 = sumstats$a0,
              a1 = sumstats$a1
              ),
          1,
          paste,
          collapse = '_')

  # extract position and chromosome
  locus <- do.call(rbind, strsplit(sumstats$locus, split = "\\:"))
  colnames(locus) <- c('chr','pos')
  sumstats <- cbind(sumstats, locus)

  # calculate MAF
  sumstats$MAF <- unlist(lapply(sumstats$AF, function(af) min(af, 1-af)))
 
  if (trait %in% "binary"){
     n_case <- as.numeric(unique(sumstats$n_cases))
     n_control <- as.numeric(unique(sumstats$n) - n_case)
     stopifnot(length(c(n_case, n_control)) == 2)
     n_eff_final <- 4 / ((1 / n_control) + (1 / n_case))
  } else {
     n_eff_final <- sumstats$n
  }
  
  # re-create clean summary stat matrix 
  sumstats <-
      data.table(
          chr = sumstats$chr,
          pos = as.integer(sumstats$pos),
          rsid = sumstats$rsid,
          a1 = sumstats$a1,
          a0 = sumstats$a0,
          n_eff = n_eff_final,
          beta_se = sumstats$standard_error,
          p = sumstats$p_value,
          beta = sumstats$beta,
          INFO = NA,
          MAF = sumstats$MAF
      )
    
   
    # remove SNPs that did not converge and thus has NAs in beta values.
    convergence <- !is.na(sumstats$beta)
    n_unconverged <- sum(!convergence)
    n <- length(convergence)
    if (n_unconverged > 0 & remove_failed){
        sumstats <- sumstats[convergence,]  
        write(paste("[Summary statistics]: Removed", n_unconverged,
                    "of", n,"variants that did not converge"),stderr())
    }

    return(sumstats)
}

