#' @title evaluate genotype agreement between whatshap and shapeit5e
#' @param dt a data.table with phased sets (GT_rb) and orther phased genotypes (GT) through S5.
#' @param labels label bins that match dt$bin
#' @param pp_cutoff phasing confidence cutoff
#' @export

# a function to evaluate the phasing accuracy by lab
eval_gt_agreement_by_bin <- function(dt, labels, pp_cutoff = NULL){
    stopifnot("GT" %in% colnames(dt))
    stopifnot("GT_rb" %in% colnames(dt))
    stopifnot("GT_rb" %in% colnames(dt))
    stopifnot("PP" %in% colnames(dt)) # phasing quality from S5
    stopifnot("rsid" %in% colnames(dt))
    stopifnot("bin" %in% colnames(dt))
    by_bin <- lapply(labels, function(lab){
        bool_rsid <-  dt$bin %in% lab
        if (!is.null(pp_cutoff)){
            stopifnot(pp_cutoff <= 1)
            bool_rsid <- bool_rsid & dt$PP >= pp_cutoff
        }
        variants <- unique(dt$rsid[bool_rsid])
        phased_sets <- unique(dt$PS_rb[dt$rsid %in% variants])
        by_ps <- lapply(phased_sets, function(ps){
            bool_sample <- (dt$PS_rb %in% ps)
            samples <- unique(dt$s[bool_sample])
            by_sample <- lapply(samples, function(sid){
                bool_ps_sid <- (dt$s %in% sid) & (dt$PS_rb %in% ps)
                dt_ps_sid <- dt[bool_ps_sid,]
                if (any(variants %in% dt_ps_sid$rsid) & (nrow(dt_ps_sid) == 2)){
                    gt_rb <- dt_ps_sid$GT_rb
                    gt_rb_mirror <- mirror_gt(dt_ps_sid$GT_rb)
                    gt <- dt_ps_sid$GT
                    match <- any(all(gt == gt_rb) | all(gt == gt_rb_mirror))
                    min_mac <- min(dt_ps_sid$AC)
                    max_mac <- max(dt_ps_sid$AC)
                    min_pp <- min(dt_ps_sid$PP, na.rm = TRUE)
                    max_pp <- max(dt_ps_sid$PP, na.rm = TRUE)
                    result <- data.frame(
                        sid = sid,
                        ps = ps,
                        match = match,
                        min_mac = min_mac,
                        max_mac = max_mac,
                        min_pp = min_pp,
                        max_pp = max_pp,
                        bin = lab
                    )
                    return(result)
                } 
            })  
            by_sample <- do.call(rbind, by_sample)
            return(by_sample)
        })
        by_ps <- do.call(rbind, by_ps)
        return(by_ps)
    })
    final <- do.call(rbind, by_bin)
    return(final)
}



