
library(argparse)
library(data.table)

main <- function(args){

    spa <- fread(args$in_spa_file)
    positions <- fread(args$coordinates)
    mrg <- merge(spa, positions, all.x = TRUE, by.x = 'MarkerID', by.y = 'ensembl_gene_id')
        
    mrg$fdr <- stats::p.adjust(mrg$p.value, method = 'fdr')
    mrg <- mrg[order(mrg$fdr),]
    mrg <- mrg[mrg$fdr < args$fdr_cutoff,]

    out <- data.table(
        gene = gsub("\\..+$",'',mrg$ensembl_gene_id),
        fdr = mrg$fdr,
        contig = mrg$chromosome_name,
        start = mrg$start_position - as.numeric(args$flanking_bp),
        end = mrg$end_position + as.numeric(args$flanking_bp)
     )

    if (nrow(out) > 0){
        outfile = paste0(args$out_prefix, '.tsv.gz')
        write(paste("writing to", outfile), stderr())
        fwrite(out, outfile, sep = '\t')
    } else {
        write(paste("No genes signifcant at FDR <",args$fdr_cutoff),stderr())
    }
}

parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, help = "Where should the results be written?")
parser$add_argument("--in_spa_file", default=NULL, help = "in directory")
parser$add_argument("--coordinates", default=NULL, help = "prefix of files to be loaded")
parser$add_argument("--flanking_bp", default=0, help = "prefix of files to be loaded")
parser$add_argument("--fdr_cutoff", default=0.20, help = "Path to list of phenotypes to be loaded.")
args <- parser$parse_args()

main(args)


