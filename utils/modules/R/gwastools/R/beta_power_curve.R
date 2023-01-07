#' @title calculate ranges for detection
#' @param cases integer with number of cases in trait
#' @param controls integer with number of controls
#' @param length.out how finegrained should the curve be?
#' @param pw.tresh power threshold
#' @param p.treshold p-value threshold
#' @param maf_lower lower bound for minor allele frequency
#' @param maf_upper upper bound for MAF
#' @param beta_lower lower bound for effect size
#' @param beta_upper upper bound for effect size
#' @return a vector of inverse rank noramlized phenotype values.
#
# based on https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS3.html


beta_power_curve <- function(cases, controls, length.out = 500,
                             pw.thresh = 0.8, p.threshold = 0.05, 
                             maf_lower = 3/175000, maf_upper = 0.05,
                             beta_lower = 0, beta_upper = 100){
    
    
    stopifnot(maf_lower < maf_upper)
    stopifnot(beta_lower < beta_upper)
    q = qchisq(p.threshold, df = 1, lower = F) #chi-square value corresp. significance threshold
    f = seq(maf_lower, maf_upper, length = length.out)
    b = seq(beta_lower, beta_upper, length = length.out)
    
    pw = rep(NA, length(b)) #power at each candidate b
    b.for.f = rep(NA,length(f)) #for each f gives the b value that leads to target power
    for(i in 1:length(f)){ 
        ncp = cases*controls / (cases+controls)*2*f[i]*(1-f[i])*b^2
        pw = pchisq(q, df = 1, ncp = ncp, lower = F)
        b.for.f[i] = suppressWarnings(b[ min( which(pw > pw.thresh) ) ])
    }
    return(data.table(f=f, beta.for.f=b.for.f))
}


