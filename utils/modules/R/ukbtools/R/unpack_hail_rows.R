# a function to extract hail rows (specifically INFO and VEP field)
unpack_hail_rows <- function(data){
    info <- format_from_hail(dfv$info)
    vep <- format_from_hail(dfv$vep)
    #dbnsfp <- format_from_hail(dfv$dbnsfp)
    mat <- setDT(cbind(
        dfv[,c('locus','alleles','snpid','rsid')],
        extract_hail_field(info),
        extract_hail_field(vep)
        #extract_hail_field(dbnsfp)
    ))
    return(mat)
}

