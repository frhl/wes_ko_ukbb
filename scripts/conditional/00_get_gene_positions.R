## NOTE: THIS CODE HAS BEEN RUN LOCALLY AND THE RESULTING FILE UPLOADED
library(biomaRt)
library(data.table)
biolist <- as.data.frame(listMarts())
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
t2g <- getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = ensembl)
t2g <- t2g[t2g$chromosome_name %in% c(1:22,'X','Y'),]
fwrite(t2g, file = '211124_ensgid_to_grch38_pos.tsv.gz', quote = FALSE, sep = '\t', row.names = FALSE)


