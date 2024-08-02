#!/usr/bin/env Rscript
library(data.table)
library(argparse)
library(stringr)
library(digest)

# function to extract annotation, af, pp and mode
str_extract_annotation <- function(x){
    return(stringr::str_extract(x, pattern=("(pLoF_damaging_missense)|(pLoF)|(damaging_missense)|(synonymous)|(non_coding)")))
}

# extract AF cutoff
str_extract_af <- function(x){
    return(as.numeric(gsub("af","",stringr::str_extract(x, pattern=("af(0\\.)?[0-9]+"))))/100)
}

# extract mode
str_extract_mode <- function(x){
    return(stringr::str_extract(x, pattern=("(recessive)|(additive)")))
}

# extract pp-value
str_extract_pp <- function(x){
  return(as.numeric(gsub("pp","",stringr::str_extract(x, pattern=("pp(0\\.)?[0-9]+")))))
}

# get index of basename
str_extract_basename_index <- function(x, index, delim="\\."){
  return(unlist(strsplit(basename(x), split=delim))[index])
}

# extract mode
str_extract_sample_subset <- function(x){
    return(stringr::str_extract(x, pattern=("(keep.wes200k)|(remove.wes200k)")))
}



main <- function(args){

    print(args)

    input_dir <- args$input_dir
    gene_map_path <- args$gene_map_path
    pattern <- args$pattern
    out_prefix <- args$out_prefix

    files <- list.files(input_dir, full.names=TRUE, pattern=pattern)
    stopifnot(length(files)>0)
    n_files <- length(files)
    write(paste0("Reading ", n_files," files..."), stderr())
    print(head(files))

    fread_saige_alt <- function(file){
        write(paste("Reading", basename(file)), stderr())
        dt <- fread(file)
        dt$pp_cutoff <- str_extract_pp(file)
        dt$af_cutoff <- str_extract_af(file)
        dt$annotation <- str_extract_annotation(file)
        dt$mode <- str_extract_mode(file)
        dt$trait <- str_extract_basename_index(file, 3)
        dt$fname_hash <- substr(digest(file, algo = "sha1"),start=1, stop=6)
        dt$sample_subset <- str_extract_sample_subset(file)
        return(dt)
    }
    
    # merge
    write("Merging into single data.table", stderr())
    lst <- lapply(files, fread_saige_alt)
    merged <- do.call(rbind, lst)
    
    # check if chrom in contig list
    if (!all(grepl(merged$CHR, pattern="chr"))){
        merged$CHR <- paste0("chr", merged$CHR)
    }    

    # get hgnc symbol
    if (!is.null(gene_map_path)){
        map_df <- fread(gene_map_path)
        ensgid_to_hgnc <- map_df$hgnc_symbol
        names(ensgid_to_hgnc) <- map_df$ensembl_gene_id
        merged$hgnc_symbol <- ensgid_to_hgnc[merged$MarkerID]
    }

    # 
    out_file <- paste0(out_prefix, ".txt.gz")
    fwrite(merged, out_file, sep = "\t")
}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_dir", default=NULL, required = TRUE, help = "path to input directory")
parser$add_argument("--pattern", default=NULL, required = TRUE, help = "pattern for input")
parser$add_argument("--gene_map_path", default=NULL, required = FALSE, help = "path gene map, e.g. 20524_hgnc_ensg_enst_chr_pos.txt.gz")
parser$add_argument("--column_info", default=NULL, required = FALSE, help = "Add a column with content")
parser$add_argument("--column_info_names", default=NULL, required = FALSE, help = "Names to be assigned")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

