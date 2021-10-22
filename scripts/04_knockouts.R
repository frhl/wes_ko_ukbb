# module purge
# conda activate rpy
# Rscript 02_knockouts.R --in_dir ../derived/knockouts/all/211013_ptv

# setup paths and libs
library(argparse)
library(data.table)
library(ggplot2)
library(dplyr)

#setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, help = "What directory should files be loaded from?")
parser$add_argument("--in_pattern", default = 'knockouts', help = "what string should the file contains?")
parser$add_argument("--print_input", default= FALSE, help = "print the input files headers to the terminal")
args <- parser$parse_args()

# read in .bgz files 
zcat_fread <- function(path,...){
	cmd = paste("zcat", path)
	return(fread(cmd = cmd,...))
}

# get directory to protein coding genes
protein_coding <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/protein_coding_genes.tsv')
protein_coding <- protein_coding$ensembl_gene_id[protein_coding$gene_biotype == 'protein_coding']

# get knockouts
ko_files = sort(list.files(args$in_dir, full.names = TRUE, pattern = args$in_pattern))
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
#write("## Done ##",stdout())

total_samples=175000
total_genes=length(unique(protein_coding))
# setup overall stats
samples = unique(dt$s)
n_samples = length(samples)
genes = unique(dt$gene_id)
n_genes = length(genes)
#print(paste(n_samples,'samples and',n_genes,'protein coding genes were loaded.'))

# count categories
knockouts = length(unique(dt[dt$knockout == 1]$s))
homozygous = length(unique(dt[dt$knockout == 1 & dt$csqs == 'HO']$s))
compound_heterozygous = length(unique(dt[dt$knockout == 1 & dt$csqs == 'CH']$s))
both = length(unique(dt[dt$knockout == 1 & dt$csqs == 'CH+HO']$s))
n = length(unique(dt$s))
col_all = c(n, homozygous, compound_heterozygous, knockouts)

# knockout stats by sample
samples_pct_ho_ko = round(100*(homozygous / total_samples), 2)
samples_pct_ch_ko = round(100*(compound_heterozygous / total_samples), 2)
samples_pct_ch_ho_ko = round(100*(both / total_samples), 2)
write("\n### Knockouts by samples ###",stdout())
write(paste0(homozygous,"/",total_samples, ' (',samples_pct_ho_ko,'%) of unqiue samples are homozygous KOs'), stdout())
write(paste0(compound_heterozygous,"/",total_samples, " (",samples_pct_ch_ko,'%) of unique samples are compound heterozygous KOs'),stdout())
write(paste0(both,"/", total_samples ," (", samples_pct_ch_ho_ko,'%) of unique samples are homozygous AND compound heterozygous KOs'),stdout())
write(paste0(knockouts,"/", total_samples ," (", samples_pct_ch_ho_ko,'%) of unique samples are KOs'),stdout())

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


# what genes are affected
#n_genes = length(unique(dt$gene_id[dt$knockout==1]))
#genes = table(dt$gene_id[dt$knockout==1])
#print(paste(n_genes,'genes are found in human knockouts'))

## assessing genes involved in mendelian disorders
# load OMIM
#omim_mapping <- fread('/well/lindgren/flassen//ressources/genesets/genesets/data/omim/mim2gene.txt')
#colnames(omim_mapping) <- c('MIM','entry','entrez','hgnc_symbol','gene_id')
#omim_knockouts <- omim_mapping[omim_mapping$gene_id %in% dt$gene_id[dt$knockout==1]]

# relate to mendelian disorderes
#mim_map <- fread('/well/lindgren/flassen//ressources/genesets/genesets/data/omim/morbidmap.txt')
#colnames(mim_map) <- c('phenotype','gene_symbols','MIM','cyto_location')
#omim_merge <- merge(mim_map, omim_knockouts)

# count how many knockout genes are associated with a phenotype
#print(paste(length(unique(omim_merge$gene_id)), 'genes associated with',length(unique(omim_merge$phenotype)),'phenotypes (OMIM)'))
#print(paste(length(unique(omim_merge$gene_id)), 'genes associated with',length(unique(omim_merge$MIM)),'MIM numbers (OMIM)'))

# count how many samples are affected
#ukbb_all_in_merge <- length(unique(dt$s[dt$knockout == 1 & dt$gene_id %in% omim_merge$gene_id]))
#print(paste(ukbb_all_in_merge, 'WES samples (all) who have knockouts of OMIN genes'))

