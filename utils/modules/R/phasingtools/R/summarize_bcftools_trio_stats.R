#' @title summarize bcftools trio stats output
#' @param d a data.frame (typically file generated with bcftools --trio-stats)
#' @param alpha confidence level for binary CIs
#' @param decimals how many decimals should SERs be rounded with?
#' @export 

summarize_bcftools_trio_stats <- function(d, alpha = 0.05, decimals = 3){
    
    stopifnot(nrow(d) > 50)
    stopifnot(ncol(d) == 8)
    colnames(d) <- 
        c('TRIO',
          'PID',
          'MID',
          'IID',
          'tested',
          'errors_mendel', 
          'errors_switch',
          'errors_switch_pct')
    trios <- nrow(d)
    # get median numbers
    median_tested <- median(d$tested)
    median_errors_switch <- median(d$errors_switch)
    median_errors_mendel <- median(d$errors_mendel)
    # get sums across all trios
    tested <- sum(d$tested)
    errors_switch <- sum(d$errors_switch)
    errors_mendel <- sum(d$errors_mendel)
    # estimate binary CI
    ci <- Hmisc::binconf(errors_switch, tested, alpha = alpha)
    ci_ser_est_pct <- paste0(round(ci[1]*100, decimals),'%')
    ci_ser_error_pct <- paste0(round(abs(ci[1] - ci[3])*100, decimals),'%')
    
    return(
        data.frame(
            trios, 
            n_tested = tested,
            n_switch = errors_switch,
            n_mendel = errors_mendel,
            median_tested, 
            median_errors_switch, 
            median_errors_mendel, 
            ci_ser_est = ci[1] ,
            ci_ser_error = abs(ci[1] - ci[3]),
            ci_ser_est_pct, 
            ci_ser_error_pct
        )
    )
}




