
library(data.table)
library(argparse)
library(pscl)

main <- function(args){
   
    # load hgncid to ensembl mapping ---------
    bridge <- fread(args$path_bridge)
    bridge <- bridge[bridge$hgnc_symbol != "",]
    hgncid_to_ensembl <- bridge$ensembl_gene_id
    names(hgncid_to_ensembl) <- bridge$hgnc_id
    hgnc_to_ensembl <- bridge$ensembl_gene_id
    names(hgnc_to_ensembl) <- bridge$hgnc_symbol

    # load genCC columns and subset ----------
    gencc <- fread(args$path_gencc)
    cols_to_keep <- c("uuid","gene_curie","gene_symbol", "disease_title", "disease_original_title",
                      "classification_title", "moi_curie", "moi_title", "disease_original_curie")
    gencc <- gencc[,colnames(gencc) %in% cols_to_keep,,with=FALSE]
    gencc$ensembl_gene_id <- hgncid_to_ensembl[gencc$gene_curie]
    colnames(gencc)[1:9] <- paste0("gencc.",colnames(gencc)[1:9])
    gencc$hgnc_symbol <- gencc$gencc.gene_symbol
    gencc$gencc.gene_symbol <- NULL

    # write outfile
    outfile = paste0(args$out_prefix, ".gencc.txt.gz")
    write(paste0("writing to ", outfile), stderr())
    fwrite(gencc, outfile, sep = '\t')

    #### read omim ---------
    lines <- readLines(args$path_omim)
    lines <- lines[!grepl("^\\#", lines)]
    #genes_to_find <- unique(na.omit(d$hgnc_symbol))
    genes_to_find <- unique(bridge$hgnc_symbol)
    omim_entry <- rbindlist(lapply(genes_to_find, function(g) {
        line <- lines[grepl(g, lines)]
        if (length(line) > 0){
            # create lines of all matches with gene name
            line_list <- strsplit(line, split = "\t")
            omim_lines <- data.table(do.call(rbind, line_list))
            omim_lines$gene <- g
            all_genes <- strsplit(gsub(" ", "",omim_lines$V2), split = ",")
            # ensure that full gene name appear in line
            lines_to_keep <- unlist(lapply(all_genes, function(line) g %in% line))
            omim_lines <- omim_lines[lines_to_keep,]                     
            if (any(lines_to_keep)){
                return(omim_lines)
            } else {
                return(NULL)
            }
        }
    }))

    # prettify and ready to combine
    colnames(omim_entry) <- c("phenotype", "gene_symbols", "mim", "cyto_loc", "hgnc_symbol")
    omim_entry$cyto_loc <- NULL # don't need location
    omim_entry$mim <- NULL # don't need omim gene id
    genes <- unique(omim_entry$hgnc_symbol) 

    # combine omim lines to one gene per line ----------
    omim <- do.call(rbind, lapply(genes, function(g){
        g_omim <- omim_entry[omim_entry$hgnc_symbol == g,]
        # grab genes that are involved in all phenotypes
        genes_involved <- unlist(strsplit(paste(g_omim$gene_symbols, collapse = ", "), split = ", "))
        genes_involved <- sort(unique(genes_involved))
        # grab disease omim IDs
        mims <- unique(na.omit(stringr::str_extract(g_omim$phenotype, "[0-9]{4,10}+")))
        phenotype <- gsub("\\,\\ [0-9]+","",g_omim$phenotype)
        phenotype <- gsub("\\ \\([0-9]\\)", "", phenotype)
        phenotype <- unique(phenotype)
        # combine
        return(data.table(
            hgnc_symbol = g,
            ensembl_gene_id = hgnc_to_ensembl[g],
            phenotypes = paste(phenotype, collapse = "; "),
            mims = paste(mims, collapse = "; "),
           genes = paste(genes_involved, collapse = "; ")
        ))
    }))

    colnames(omim)[3:5] <- c("omim.phenotypes","omim.mim","omim.genes") 
    outfile = paste0(args$out_prefix, ".omim.txt.gz")
    write(paste0("writing to ", outfile), stderr())
    fwrite(omim, outfile, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_bridge", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_omim", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_gencc", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


