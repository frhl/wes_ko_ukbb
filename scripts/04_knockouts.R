# module purge
# conda activate rpy
# Rscript 04_knockouts.R --in_dir derived/knockouts/211111 --in_pattern maf00_50 --in_csq ptv_damaging_missense_knockouts --out_prefix derived/summary

# setup paths and libs
library(argparse)
library(data.table)
library(ggplot2)
library(dplyr)

setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, help = "What directory should files be loaded from?")
parser$add_argument("--in_pattern", default = 'knockouts', help = "what string should the file contains?")
parser$add_argument("--in_csq", default = '', help = "what consequence should we subset by")
parser$add_argument("--in_ext", default = 'tsv.bgz', help = "what file extension is required")
parser$add_argument("--print_input", default= FALSE, help = "print the input files headers to the terminal")
parser$add_argument("--out_prefix", default=NULL, help = "out prefix for files")
args <- parser$parse_args()


## helpers
# read in .bgz files 
zcat_fread <- function(path,...){
	cmd = paste("zcat", path)
	return(fread(cmd = cmd,...))
}

# process variants on each strand
clean_hail_list <- function(x, comma_sub = ';') { 
    x <- gsub('(\\[)|(\\])|(\\")','',x)
    x <- gsub('\\,',comma_sub,x)
    return(x)
}


## main
# get directory to protein coding genes
protein_coding <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/protein_coding_genes.tsv')
protein_coding <- protein_coding$ensembl_gene_id[protein_coding$gene_biotype == 'protein_coding']

# get knockouts
ko_files = sort(list.files(args$in_dir, full.names = TRUE, pattern = args$in_pattern))
ko_files = ko_files[grepl(args$in_csq, ko_files)]
regex = paste0('*.',args$in_ext,'$')
ko_files = ko_files[grepl(regex, ko_files)]

if (length(ko_files) < 1) stop('No files were found in the specified directory (with the specified pattern)!')


# combine knockouts in single file
write(paste("## Loading files from", args$in_dir,"##"),stdout())
dt <- setDT(do.call(rbind, lapply(ko_files, function(f){
    d = zcat_fread(f)
    d = d[d$gene_id %in% protein_coding]
    nrows = nrow(d)
    name = basename(f)
    chr = as.numeric(gsub('chr','',unlist(lapply(strsplit(name, split = '_'), function(x) x[7]))))
    d$chr = chr
    carriers = length(unique(d$s))
    kos = length(unique(d$s[d$knockout == 1]))
    ko_genes = length(unique(d$gene_id[d$knockout == 1]))
    out = paste0("loaded ",name," summary={lines = ",nrows,"; carriers = ",carriers,"; KO samples = ",kos,"; KO genes = ",ko_genes,"}")
    write(out,stdout())
    # print out head
    dp = d[!is.na(d$phase1) & !is.na(d$phase2),]
    if (args$print_input) print(head(dp))
    return(d)
})))

# Summarize results
total_samples=176929
total_genes=length(unique(protein_coding))
# setup overall stats
samples = unique(dt$s)
n_samples = length(samples)
genes = unique(dt$gene_id)
n_genes = length(genes)

# count categories
knockouts = length(unique(dt[dt$knockout == 1]$s))
homozygous = length(unique(dt[dt$knockout == 1 & dt$csqs == 'HO']$s))
compound_heterozygous = length(unique(dt[dt$knockout == 1 & dt$csqs == 'CH']$s))
both = length(unique(dt[dt$knockout == 1 & dt$csqs == 'CH+HO']$s))
n = length(unique(dt$s))
col_all = c(n, homozygous, compound_heterozygous, knockouts)

# count categories per individual
sample_ko = dt[dt$knockout == 1]$s
ko_per_individual = length(sample_ko) / length(unique(dt$s))
print(ko_per_individual)

# knockout stats by sample
samples_pct_ho_ko = round(100*(homozygous / total_samples), 2)
samples_pct_ch_ko = round(100*(compound_heterozygous / total_samples), 2)
samples_pct_ch_ho_ko = round(100*(both / total_samples), 2)
samples_pct_ko = round(100*(knockouts / total_samples), 2)
write("\n### Knockouts by samples ###",stdout())
write(paste0(homozygous,"/",total_samples, ' (',samples_pct_ho_ko,'%) of unqiue samples are homozygous KOs'), stdout())
write(paste0(compound_heterozygous,"/",total_samples, " (",samples_pct_ch_ko,'%) of unique samples are compound heterozygous KOs'),stdout())
write(paste0(both,"/", total_samples ," (", samples_pct_ch_ho_ko,'%) of unique samples are homozygous AND compound heterozygous KOs'),stdout())
write(paste0(knockouts,"/", total_samples ," (", samples_pct_ko,'%) of unique samples are KOs'),stdout())

# knockout stats by gene
genes_knockout = length(unique(dt$gene_id[dt$knockout == 1])) 
genes_ho = length(unique(dt$gene_id[dt$csqs == 'HO'])) 
genes_ch = length(unique(dt$gene_id[dt$csqs == "CH"])) 
genes_pct_ho_ko = round(100*(genes_ho / total_genes), 2)
genes_pct_ch_ko = round(100*(genes_ch / total_genes), 2)
genes_pct_ko = round(100*(genes_knockout / total_genes), 2)

write("\n### Knockouts by Genes ###",stdout())
write(paste0(genes_ho,"/",total_genes, ' (',genes_pct_ho_ko,'%) of unique genes are involved in homozygous KOs'), stdout())
write(paste0(genes_ch,"/",total_genes, " (",genes_pct_ch_ko,'%) of unique genes are involved in compound heterozygous KOs'),stdout())
write(paste0(genes_knockout,"/", total_genes ," (", genes_pct_ko,'%) of unique genes are involed in KOs'),stdout())

# Write tables of aggregated variants
dt$phase1 <- clean_hail_list(dt$phase1)
dt$phase2 <- clean_hail_list(dt$phase2)

# if phase is NA replace with "."
dt$phase1[is.na(dt$phase1)] <- '.'
dt$phase2[is.na(dt$phase2)] <- '.'

# check if involved in OMIM
omim <- fread('/well/lindgren/flassen//ressources/genesets/genesets/data/omim/211103_morbidmap_by_gene.txt')

# map ensgid to hgnc_symbol
hgnc_link <- fread('/well/lindgren/flassen//ressources/genesets/genesets/data/hgnc/211026_hgnc_ensgid_link.csv')
dt$in_omim <- dt$gene_id %in% omim$ensgid
dt <- merge(dt, hgnc_link, by.x = 'gene_id', by.y = 'ensgid', all.x = TRUE)
dt <- cbind(dt[,c('s','gene_id','hgnc_symbol')], dt[,-c('s','gene_id','hgnc_symbol')])

# write out table
outfile = paste0(args$out_prefix, '.tsv.gz')
fwrite(dt, outfile, sep = '\t', row.names = FALSE)

# Write out unique genes involved in OMIM
genes_omim <- omim$ensgid
genes_knockout_omim = length(unique(dt$gene_id[dt$knockout == 1 & dt$gene_id %in% genes_omim])) 
genes_ho_omim = length(unique(dt$gene_id[dt$csqs == 'HO' & dt$gene_id %in% genes_omim])) 
genes_ch_omim = length(unique(dt$gene_id[dt$csqs == "CH" & dt$gene_id %in% genes_omim])) 

genes_omim_pct_ho_ko = round(100*(genes_ho_omim / total_genes), 2)
genes_omim_pct_ch_ko = round(100*(genes_ch_omim / total_genes), 2)
genes_omim_pct_ko = round(100*(genes_knockout_omim / total_genes), 2)

write("\n### Knockouts by OMIM Genes ###",stdout())
write(paste0(genes_ho_omim,"/",total_genes, ' (',genes_omim_pct_ho_ko,'%) of unique genes are involved in homozygous OMIM KOs'), stdout())
write(paste0(genes_ch_omim,"/",total_genes, " (",genes_omim_pct_ch_ko,'%) of unique genes are involved in compound heterozygous OMIM KOs'),stdout())
write(paste0(genes_knockout_omim,"/", total_genes ," (", genes_omim_pct_ko,'%) of unique genes are involed in OMIM KOs'),stdout())



