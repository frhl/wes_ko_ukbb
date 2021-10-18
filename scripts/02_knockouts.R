
# setup paths and libs
library(data.table)
library(ggplot2)
library(dplyr)
setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')

# protein coding genes
protein_coding <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/protein_coding_genes.tsv')
protein_coding <- protein_coding$ensembl_gene_id[protein_coding$gene_biotype == 'protein_coding']

# get knockouts
ko_files = list.files('derived/knockouts/all/211013_ptv/', full.names = TRUE, pattern = 'knockouts')

# combine knockouts in single file
dt <- setDT(do.call(rbind, lapply(ko_files, function(f){
    d = fread(f)
    d = d[d$gene_id %in% protein_coding]
    name = basename(f)
    chr = as.numeric(gsub('chr','',unlist(lapply(strsplit(name, split = '_'), function(x) x[7]))))
    d$chr = chr
    return(d)
})))

# setup overall stats
samples = unique(dt$s)
genes = unique(dt$gene_id)
print(paste(length(samples),'samples and',length(genes),'protein coding genes'))

# count categories
knockouts = length(unique(dt[dt$knockout == 1]$s))
homozygous = length(unique(dt[dt$knockout == 1 & dt$csqs == 'HO']$s))
compound_heterozygous = length(unique(dt[dt$knockout == 1 & dt$csqs == 'CH']$s))
both = length(unique(dt[dt$knockout == 1 & dt$csqs == 'CH+HO']$s))
n = length(unique(dt$s))
col_all = c(n, homozygous, compound_heterozygous, knockouts)

# combine in data.frame
df <- data.frame(col_all)
rownames(df) <- c('n','homozygous knockouts', 'compound heterozygous knockouts', 'total knockouts')
print(df)

# what genes are affected
n_genes = length(unique(dt$gene_id[dt$knockout==1]))
genes = table(dt$gene_id[dt$knockout==1])
print(paste(n_genes,'genes are found in human knockouts'))

## assessing genes involved in mendelian disorders
# load OMIM
omim_mapping <- fread('/well/lindgren/flassen//ressources/genesets/genesets/data/omim/mim2gene.txt')
colnames(omim_mapping) <- c('MIM','entry','entrez','hgnc_symbol','gene_id')
omim_knockouts <- omim_mapping[omim_mapping$gene_id %in% dt$gene_id[dt$knockout==1]]

# relate to mendelian disorderes
mim_map <- fread('/well/lindgren/flassen//ressources/genesets/genesets/data/omim/morbidmap.txt')
colnames(mim_map) <- c('phenotype','gene_symbols','MIM','cyto_location')
omim_merge <- merge(mim_map, omim_knockouts)

# count how many knockout genes are associated with a phenotype
print(paste(length(unique(omim_merge$gene_id)), 'genes associated with',length(unique(omim_merge$phenotype)),'phenotypes (OMIM)'))
print(paste(length(unique(omim_merge$gene_id)), 'genes associated with',length(unique(omim_merge$MIM)),'MIM numbers (OMIM)'))

# count how many samples are affected
ukbb_all_in_merge <- length(unique(dt$s[dt$knockout == 1 & dt$gene_id %in% omim_merge$gene_id]))
print(paste(ukbb_all_in_merge, 'WES samples (all) who have knockouts of OMIN genes'))



