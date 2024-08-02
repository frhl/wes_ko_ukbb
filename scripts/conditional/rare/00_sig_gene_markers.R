
# Takes gene hits that are significant in primary analysis
# and digs out all variants in the gene as annotated by VEP.

library(argparse)
library(data.table)

main <- function(args){


    if (!file.exists(args$in_spa_file)) stop(paste(args$in_spa_file, "does not exists!"))
    stopifnot(grepl("CHR", args$vep_path_with_CHR))
    spa <- fread(args$in_spa_file)
    spa_n <- nrow(spa)

    # test P-value
    p_cutoff <- as.numeric(args$p_cutoff)
    stopifnot(p_cutoff > 0 && p_cutoff < 1)
    spa <- spa[spa$p.value < p_cutoff,]
    write("Note: using standard column $p.value to check cutoff", stderr())

    vep_path_old <- ""
    vep <- NULL

    # only continue if something is significant
    if (nrow(spa) > 0){
        
        # iterate over each row and recover variants for the gene
        lst <- lapply(1:nrow(spa), function(row_idx){
        
            row <- spa[row_idx,]
            gene <- row$MarkerID
            chrom <- row$CHR

            # read in vep
            vep_path <- gsub("chrCHR", chrom, args$vep_path_with_CHR)
            stopifnot(file.exists(vep_path))
            if (vep_path != vep_path_old){
                vep <<- fread(vep_path)
            }

            # save old VEP path as global
            vep_path_old <<- vep_path
            
            # filter to genes in that are significant
            vep_subset <- vep[vep$csqs.gene_id %in% gene,]
            vep_subset <- vep_subset[,c("varid","csqs.gene_id","consequence_category")]
            colnames(vep_subset) <- c('rsid',"ensembl_gene_id","consequence_category")
            return(vep_subset)
        
        })
        
        # combine variants for phenotype
        d <- do.call(rbind, lst)
        d$phenotype <- args$phenotype
        stopifnot(nrow(d) > 0)
        
        # write file
        outfile = paste0(args$out_prefix, '.txt.gz')
        write(paste("writing to", outfile), stdout())
        fwrite(d, outfile, sep = '\t')
        
    } else {
        write(paste0("No lines in SPA file after subsetting by P-value for ", args$in_spa_file), stderr())
    }

}

parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, help = "Where should the results be written?")
parser$add_argument("--in_spa_file", default=NULL, help = "In file to be assessed")
parser$add_argument("--vep_path_with_CHR", default=NULL, help = "Path to VEP file with 'chrCHR' format so that it can be gsubbed")
parser$add_argument("--p_cutoff", default=5e-6, help = "Path to list of phenotypes to be loaded.")
parser$add_argument("--phenotype", default=NULL, help = "Current phenotype")
args <- parser$parse_args()

main(args)


