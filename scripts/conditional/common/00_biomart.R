setwd('/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb')

library(biomaRt)
library(data.table)
biolist <- as.data.frame(listMarts())
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh = NULL) # GRCh38
d <- getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position', 'strand'), mart = ensembl)
d <- d[d$chromosome_name %in% c(1:22,'X','Y'),]
d$upstream <- factor(ifelse(d$strand == 1, 'LOWER', 'HIGHER'))
outfile = 'data/genes/220310_ensgid_grch38_pos.tsv.gz'
fwrite(d, outfile, quote = FALSE, sep = '\t', row.names = FALSE)


