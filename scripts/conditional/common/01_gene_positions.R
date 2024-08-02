
library(argparse)
library(data.table)

main <- function(args){

    stopifnot(file.exists(args$in_spa_file))
    stopifnot(file.exists(args$coordinates))

    spa <- fread(args$in_spa_file)
    spa_n <- nrow(spa)
    positions <- fread(args$coordinates)
    mrg <- merge(spa, positions, all.x = TRUE, by.x = 'MarkerID', by.y = 'ensembl_gene_id')

    out <- data.table(
        gene = gsub("\\..+$",'',mrg$ensembl_gene_id),
        pvalue = mrg$p.value,
        contig = mrg$chromosome_name,
        start = mrg$start_position - as.numeric(args$flanking_bp),
        end = mrg$end_position + as.numeric(args$flanking_bp),
        prs = grepl("locoprs", args$in_spa_fil)
     )

    # filter to signifcant entries
    p_cutoff <- as.numeric(args$p_cutoff)
    out <- out[out$pvalue < p_cutoff,]
    n <- nrow(out)
    msg <- paste("Subsetted", args$phenotype, "to", n, "significant entries in", args$in_spa_file)
    write(msg, stdout())


    if (n > 0){
        outfile = paste0(args$out_prefix, '.tsv.gz')
        write(paste("writing to", outfile), stderr())
        fwrite(out, outfile, sep = '\t')
    } 
}

parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, help = "Where should the results be written?")
parser$add_argument("--in_spa_file", default=NULL, help = "in directory")
parser$add_argument("--coordinates", default=NULL, help = "prefix of files to be loaded")
parser$add_argument("--flanking_bp", default=0, help = "prefix of files to be loaded")
parser$add_argument("--fdr_cutoff", default=0.20, help = "Path to list of phenotypes to be loaded.")
parser$add_argument("--p_cutoff", default=5e-6, help = "Path to list of phenotypes to be loaded.")
parser$add_argument("--phenotype", default=NULL, help = "Current phenotype")
args <- parser$parse_args()

main(args)

