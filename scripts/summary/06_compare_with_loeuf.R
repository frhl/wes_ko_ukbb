# module purge
# conda activate rpy
# Rscript  06_compare_with_loeuf.R  --out_prefix derived/plots/ukb_wes200k_knockout_loeuf --MAF 00_01 --CATEGORY ptv_damaging_missense

library(argparse)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggrastr)

# setup paths and load helpers
setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')
devtools::load_all('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/ukbtools')

parser <- ArgumentParser()
parser$add_argument("--out_prefix", default='derived/plots/ukb_wes200k_knockout_loeuf', help = "Dest dir")
parser$add_argument("--MAF", default='00_01', help = "Minor allele threshold in format XX_YY")
parser$add_argument("--CATEGORY", default='ptv_damaging_missense', help = "consequence cateogories")
args <- parser$parse_args()

# setup parameters
#MAF = '00_01'
#CATEGORY = 'ptv_damaing_missense'
MAF = args$MAF
CATEGORY = args$CATEGORY
TITLE=paste0(MAF,'_',CATEGORY)
print(TITLE)

# Load data
knockouts <- list.files('derived/knockouts/211111/', full.names = TRUE)
dt <- load_knockout_bundle(knockouts, MAF, CATEGORY)
dt <- melt(dt, id.vars = 'gene_id')
colnames(dt) <- c('gene','mutation','samples')
dt <- dt[dt$samples > 0]

# translate genes into hgnc_symbol
hgnc_bridge <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/hgnc/211026_hgnc_ensgid_link.csv')
colnames(hgnc_bridge)[2] <- 'gene'
dt <- merge(dt, hgnc_bridge, by = 'gene', all.x = TRUE)

# Load gnomAD LOEUF data
constraints <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv')
constraints$loeuf <- constraints$oe_lof_upper 
constraints <- constraints[constraints$canonical == TRUE,]

# create deciles
deciles_loeuf_seq <- seq(0,1,by = 0.1)
deciles_loeuf <- quantile(constraints$loeuf, probs = deciles_loeuf_seq, na.rm = T)
constraints$decile <- cut(constraints$loeuf, deciles_loeuf)
levels(constraints$decile) <- deciles_loeuf_seq*100
constraints$decile_plot <-  paste0(constraints$decile,'-', as.numeric(as.character(constraints$decile))+10)
constraints$decile_plot[constraints$decile_plot == 'NA-NA'] <- NA

# merge with LOEUF
dt <- merge(dt, constraints[,c('gene','decile_plot','pLI','oe_lof_upper', 'oe_lof_lower')], by.x = 'hgnc_symbol', by.y = 'gene', all.x = TRUE)

# remove heterozygous variants
dt <- dt[dt$mutation != 'HE',]

# remove genes that is knocked out in more than 50% of samples
genes <- dt$gene[dt$samples > 100000]
print(dt[dt$gene %in% genes,])
print(paste('Note: excluded',length(genes),'knockout genes in more than 50% of samples:'))
print(paste(genes, collapse = ','))
dt <- dt[! dt$gene %in% genes]

# setup decile histogram
plt_hist <- ggplot(dt, aes(x=decile_plot, y=samples, fill=mutation, label=gene)) +
    geom_bar(stat='identity') + 
    xlab('LOEUF Decile') +
    ylab('Samples (n)') +
    ggtitle('knockouts by LOEUF decile', TITLE)

# label genes with most samples in each group
genes_ok <- na.omit(unlist(lapply(unique(dt$decile_plot), function(x){
    d <- dt[dt$decile_plot %in% x,]
    d <- d[rev(order(d$samples)),]
    d <- head(d$hgnc_symbol, 3)
    return(d)
})))

dt$label <- ''
dt$label[dt$hgnc_symbol %in% genes_ok] <-  dt$hgnc_symbol[dt$hgnc_symbol %in% genes_ok]

# generate jitter plot with top 3 genes in each group labelled
pos <- position_jitter(width = 0.3, seed = 2)
plt_jitter <- ggplot(dt[!is.na(dt$decile_plot),], aes(x=decile_plot, y=samples, color=mutation, label=label)) +
    geom_jitter_rast(position = pos, size = 2) +     
    geom_text_repel(
        color = 'black', 
        position = pos, 
        hjust = "outward", 
        direction = "y",
        max.overlaps = 200,
        min.segment.length = unit(0.25, 'lines')
    ) + 
    xlab('LOEUF Decile') +
    ylab('Samples (n)') +
    ggtitle('Knockouts by LOEUF deciles', TITLE)


# save jitter
outfile <- paste0(args$out_prefix, '_jitter_', TITLE) 
write(paste('saving to', outfile),stdout())
ggsave(filename = paste0(outfile,'.pdf'), plt_jitter, height = 8, width = 17)
fwrite(dt[!is.na(dt$decile_plot),], paste0(outfile,'.tsv'), row.names = FALSE, sep = '\t')

# save hist
outfile <- paste0(args$out_prefix, '_hist_', TITLE) 
write(paste('saving to', outfile),stdout())
ggsave(filename = paste0(outfile,'.pdf'), plt_hist, height = 8, width = 17)
fwrite(dt[!is.na(dt$decile_plot),], paste0(outfile,'.tsv'), row.names = FALSE, sep = '\t')




