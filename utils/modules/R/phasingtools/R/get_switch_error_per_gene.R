#' @title get switch errors per gene
#' @param switches vector of switches
#' @param pos vector of positions
#' @param chrom chromosome in the format "21" "22" etc.
#' @param gene_path path to biomat gene files.
#' @export

get_switch_error_per_gene <- function(pos, switches, chrom, gene_path = "/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/211124_ensgid_to_grch38_pos.tsv.gz"){
    
    stopifnot(is.numeric(pos))
    stopifnot(switches %in% c(0,1))
    
    # leverage genes downloaded form biomart GRCh38
    genes <- fread(gene_path)
    print(head(genes))
    stopifnot('chromosome_name' %in% colnames(genes))
    stopifnot('start_position' %in% colnames(genes))
    stopifnot("end_position" %in% colnames(genes))
    genes <- genes[genes$chromosome_name %in% as.character(chrom)]
    
    if (nrow(genes) < 1) stop(paste(chrom, "is not a valid chromosome!"))
    
    # iterate over sites
    n <- nrow(genes)
    lst <- lapply(1:n, function(idx){
        gene <- genes[idx, ]
        start <- gene$start_position
        end <- gene$end_position
        sites_boolean <- pos >= start & pos <= end
        sites_in_gene <- sum(sites_boolean)
        erors_boolean <- sites_boolean & switches > 0
        errors_in_gene <- sum(erors_boolean)
        out <- data.frame(
            gene=gene$ensembl_gene_id,
            sites=sites_in_gene,
            errors=errors_in_gene
        )
        return(out)
    })
    
    phased_genes <- do.call(rbind, lst)
    return(phased_genes)
}


