# extract a field and combine it by gene x individual in a single cell                                
extract_field_by_gene_sample <- function(df, field = 'most_severe_consequence', delim = ';'){
    samples <- unique(df$sample)
    genes <- unique(df$gene)
    stopifnot(field %in% colnames(df))
    as.data.table(do.call(rbind, lapply(samples, function(s) 
        do.call(rbind, lapply(genes, function(g){
            outdf <- df[df$sample %in% s & df$gene %in% g, ]
            outdf <- data.frame(table(outdf[[field]]))
            if (nrow(outdf) > 0){
                return(data.table(s = s, gene = g, entry = paste0(outdf$Var1, collapse = delim)))
                }
    })))))
}

