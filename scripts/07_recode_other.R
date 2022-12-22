#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(stringr)

# create mapping to used by map_by_variant
create_mapping <- function(dt, what, id="id"){
    stopifnot(what %in% colnames(dt))
    stopifnot(id %in% colnames(dt))
    mapping <- dt[[what]]
    names(mapping) <- dt[[id]]
    return(mapping)
}

# create ID based on variant id and gene id
create_id <- function(dt){
    stopifnot("gene_id" %in% colnames(dt))
    stopifnot("varid" %in% colnames(dt))
    n <- nrow(dt)
    ids <- unlist(lapply(1:n, function(idx){
        variant <- dt$varid[idx]
        gene <- dt$gene_id[idx]
        at_least_two <- stringr::str_detect(variant, ";")
         if (at_least_two){
            id <- unlist(strsplit(variant, ";"))
            id <- paste(gene, id, sep = ":")
            return(paste0(id, collapse = ";"))
        } else {
            id <- paste(gene, variant, sep = ":")
            return(id)
        }
    }))
    return(ids)
}

# map by variant
map_by_variant <- function(ids, mapping, allow_dups=TRUE, message=""){
    stopifnot(any(names(mapping) %in% ids))
    if (nchar(message)>0) write(message,stdout())
    mapped <- unlist(lapply(ids, function(v){
        at_least_two <- stringr::str_detect(v, ";")
        if (at_least_two){
            vs <- unlist(strsplit(v, ";"))
            vs <- as.character(mapping[vs])
            vs <- na.omit(vs)
            if (!allow_dups) vs <- unique(vs)
            if (length(vs)==0) return(NA)
            return(paste0(vs, collapse = ";"))
        } else {
            return(paste0(mapping[v], collapse = ";"))
        }
    }))
    return(mapped)
}



main <- function(args){

    input_path <- args$input_path
    vep_path <- args$vep_path
    out_prefix <- args$out_prefix
    csqs <- args$csqs

    dt <- fread(input_path)
    annotation <- fread(vep_path)

    # create ID fields for mapping
    colnames(annotation) <- gsub("worst_csq_by_gene_canonical\\.","",colnames(annotation))
    annotation$id <- paste0(annotation$gene_id,":", annotation$varid)
    dt$id <- create_id(dt) 

    # check that all variants are contained
    variants_in_dt <- unlist(strsplit(dt$id, split = ";"))
    ok <- as.logical(sum(variants_in_dt %in% annotation$id)/length(variants_in_dt))
    stopifnot(ok)
    
    # create mapping to VEP annotations
    map_enstid <- create_mapping(annotation, "transcript_id")
    map_exon <- create_mapping(annotation, "exon")
    map_intron <- create_mapping(annotation, "intron")
    map_revel <- create_mapping(annotation, "revel_score")
    map_cadd <- create_mapping(annotation, "cadd_phred")
    map_function_csqs <- create_mapping(annotation, "most_severe_consequence")
    map_csqs <- create_mapping(annotation, "consequence_category")
    map_codons <- create_mapping(annotation, "codons")
    map_amino_acids <- create_mapping(annotation, "amino_acids")

    # cerarte mapping to INFO annotations
    map_af <- create_mapping(annotation, "info.AF")
    map_ac <- create_mapping(annotation, "info.AC")
    map_an <- create_mapping(annotation, "info.AN")

    # map items from 
    dt$consequence_category <- csqs #" #map_by_variant(dt$id, map_csqs, message="csqs_category")
    dt$most_severe_consequence <- NA #map_by_variant(dt$id, map_function_csqs, message="most_severe_csqs")
    dt$transcript_id <- map_by_variant(dt$id, map_enstid, allow_dups = FALSE, message="trancript")
    dt$exon <- NA #map_by_variant(dt$id, map_exon, message="exon")
    dt$intron <- NA #map_by_variant(dt$id, map_intron, message="intron")
    dt$revel_score <- NA #map_by_variant(dt$id, map_revel, message="revel_score")
    dt$cadd_phred <- NA #map_by_variant(dt$id, map_cadd, message="cadd_phred")
    dt$AF <- NA #map_by_variant(dt$id, map_af, message = "af")
    dt$AC <- map_by_variant(dt$id, map_ac, message = "ac")
    dt$AN <- map_by_variant(dt$id, map_an, message = "an")
    dt$id <- NULL

    # make columns nicer
    new_order <- c("s","gene_id","transcript_id","varid","gts","AC", "AF", "AN", "hom_alt_n", "phased.a1", "phased.a2",
                   "unphased.n", "pKO" ,"knockout","consequence_category", "most_severe_consequence",
                   "revel_score", "cadd_phred",'exon','intron')
    stopifnot(all(new_order %in% colnames(dt)))
    stopifnot(all(colnames(dt) %in% new_order))
    dt <- dt[,new_order, with=FALSE]

    # create subsets and export
   
    # export all 
    outfile <- paste0(out_prefix,".",csqs,".txt.gz")
    fwrite(dt, outfile)


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "?")
parser$add_argument("--csqs", default=NULL, help = "?")
parser$add_argument("--vep_path", default=NULL, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

