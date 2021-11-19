get_knockout_sample_counts <- function(files, maf, mutation){
    
    # subset to files
    bool_maf = grepl(maf, files)
    mutation_file = paste0(mutation,'_','knockouts')
    bool_mutation = grepl(mutation_file, files)

    # load files
    files <- files[bool_mutation & bool_maf]
    stopifnot(length(files) > 0)
    combined <- setDT(do.call(rbind, lapply(files, zcat)))
    combined$csqs[combined$csqs %in% c('CH+HO','HO+CH')] <- 'HO'
    
    # combine haplotypes on string
    combined <- combined[combined$csqs %in% c('CH','HO')]
    combined$phase1 <- format_from_hail(combined$phase1)
    combined$phase2 <- format_from_hail(combined$phase2)
    combined$haplotypes <- apply(combined[,c('phase1','phase2')], 1, paste, collapse = '|')
    
    dt_ho <- combined[combined$csqs %in% c('HO')]
    dt_ch <- combined[combined$csqs %in% c('CH')]
    
    # deal with homozygous
    table_ho <- as.data.table(table(dt_ho$gene_id, dt_ho$haplotypes))
    table_ho <- table_ho[table_ho$N > 0]
    table_ho$csqs <- 'homozygous'

    # deal with compound hetz
    table_ch <- as.data.table(table(dt_ho$gene_id, dt_ho$haplotypes))
    table_ch <- table_ch[table_ch$N > 0]
    table_ch$csqs <- 'compound heterozygous'

    # combine and re-name
    table_ch_ho <- rbind(table_ho, table_ch)
    colnames(table_ch_ho) <- c('ensgid','haplotypes','sample count','csqs')

    # add chromosome and re-roder
    table_ch_ho$chr <- unlist(lapply(strsplit(table_ch_ho$haplotypes, split = '_'), function(x) gsub('chr','',x[1])))
    table_ch_ho <- table_ch_ho[order(table_ch_ho$chr),]
    return(table_ch_ho)
    
}


