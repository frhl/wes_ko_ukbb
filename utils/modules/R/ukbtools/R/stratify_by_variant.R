#s tratify by variant in long format
stratify_by_variant <- function(curdt){
   do.call(rbind, lapply(1:nrow(curdt), function(i){
    gene = curdt$gene_id[i]
    sample = curdt$s[i]
    knockout = curdt$knockout[i]
    phase1 = unlist(strsplit(format_from_hail(curdt$phase1[i]), split = ','))
    phase2 = unlist(strsplit(format_from_hail(curdt$phase2[i]), split = ','))
    outdt1 = data.table(gene = gene, sample = sample, phase = 1, knockout = knockout, snpid = phase1)
    outdt2 = data.table(gene = gene, sample = sample, phase = 2, knockout = knockout, snpid = phase2)
    return(rbind(outdt1, outdt2))
    })) 
}

