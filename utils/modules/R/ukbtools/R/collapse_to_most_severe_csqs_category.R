# collapse to
collapse_to_most_severe_csqs_category <- function(df){
    if (! 'consequence_category' %in% colnames(df)) stop('consequence_category not found!')
    snpids <- unique(df$snpid)
    do.call(rbind, lapply(snpids, function(id){
        tmp <- df[df$snpid %in% id,]
        ptv_bool <- tmp$consequence_category %in% 'ptv'
        if (sum(ptv_bool)) return(tmp[ptv_bool,])
        missense_bool <- tmp$consequence_category == 'missense'
        if (sum(missense_bool)) return(tmp[missense_bool,])
        synonymous_bool <- tmp$consequence_category == 'synonymous'
        if (sum(synonymous_bool)) return(tmp[synonymous_bool,])
        non_coding_bool <- tmp$consequence_category == 'non_coding'
        if (sum(non_coding_bool)) return(tmp[non_coding_bool,])
    }))
}


