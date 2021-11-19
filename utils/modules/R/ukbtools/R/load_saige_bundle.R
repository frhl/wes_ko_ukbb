load_saige_bundle <- function(files, maf, mutation,pheno){
    
    # set thresholds
    #maf = '00_01'
    #pheno = "Cirrhosis"
    #mutation = 'ptv_damaging_missense'

    # subset files
    pheno_mutation = paste0(mutation,'_',pheno)
    bool_maf = grepl(maf, saige_binary)
    bool_pheno = grepl(pheno, saige_binary)
    bool_mutation = grepl(pheno_mutation, saige_binary)

    # load files
    files <- saige_binary[bool_pheno & bool_mutation & bool_maf]
    stopifnot(length(files) > 0)
    combined <- setDT(do.call(rbind, lapply(files, fread)))
    
    # add analytical uniform expecation
    n <- nrow(combined)
    combined <- combined[order(combined$p.value),]
    combined$pvalue.observed <- combined$p.value
    combined$pvalue.expected <- seq(1, n)/(n + 1)
    
    # clean up
    colnames(combined)[colnames(combined)=='SNPID'] <- 'gene_id'
    
    return(combined)
    
}




