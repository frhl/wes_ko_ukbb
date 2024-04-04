library(data.table)
library(argparse)

main <- function(args){

    in_file <- args$in_file
    coordinates <- args$coordinates
    out_prefix <- args$out_prefix
    flanking_bp <- as.numeric(args$flanking_bp)
    stopifnot(file.exists(in_file))
    stopifnot(file.exists(coordinates))

    # survival file does not have correct columns yet
    d <- fread(in_file)
    d <- d[d$analysis_type == "unconditioned",]
    d$ensembl_gene_id <- d$gene
    stopifnot(nrow(d)>0)
 
    # merge file iwth positions
    positions <- fread(coordinates)
    mrg <- merge(d, positions, all.x = TRUE, by = 'ensembl_gene_id')  

    # combine out table
    out <- data.table(
        gene = gsub("\\..+$",'',mrg$ensembl_gene_id),
        pvalue = mrg$p.value,
        contig = mrg$chromosome_name,
        start = mrg$start_position - flanking_bp,
        end = mrg$end_position + flanking_bp,
        diagnosis = mrg$dianosis
    )

    # write file 
    outfile <- paste0(out_prefix, ".txt.gz")
    write(paste("writing", outfile), stdout())
    fwrite(out, outfile, sep="\t")
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_file", default=NULL, required = TRUE, help = "")
parser$add_argument("--coordinates", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
parser$add_argument("--flanking_bp", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)











