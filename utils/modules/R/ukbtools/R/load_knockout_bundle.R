load_knockout_bundle <- function(files, maf, mutation){
    
    # subset to files
    bool_maf = grepl(maf, files)
    mutation_file = paste0(mutation,'_','knockouts')
    bool_mutation = grepl(mutation_file, files)

    # load files
    files <- files[bool_mutation & bool_maf]
    stopifnot(length(files) > 0)
    combined <- setDT(do.call(rbind, lapply(files, zcat)))
 
    # aggregate knockouts type by gene id
    dt <- data.table(table(combined$gene_id, combined$csqs))
    dt <- data.table::dcast(dt, V1 ~ V2, value.var = 'N')
    colnames(dt)[1] <- 'gene_id'
    
    # clean up outout
    if ('CH+HO' %in% colnames(dt)){
        dt$HO <- dt$HO + dt$`CH+HO`
        dt$`CH+HO` <- NULL
    }
    
    # aggregate alleles by gene id
    #na_collapse <- function(x) paste0(unique(na.omit(x)), collapse = ';')
    #rsids_combined <- apply(combined[,c('phase1','phase2')], 1, na_collapse)
    #rsids_combined <- gsub('(\\[)|(\\])|(\\")', '', rsids_combined)
    #combined$rsids <- rsids_combined
    
    #dt_alleles <- combined[,c('gene_id','rsids')]
    #dt_alleles <- dt_alleles[!duplicated(dt_alleles),]

    # combine alleles and knockouts
    #dt <- merge(dt, dt_alleles, all.x = TRUE)
    

    return(dt)
}






