#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(stringr)

# create mapping to used by map_by_variant
create_mapping <- function(dt, what, id="varid"){
    stopifnot(what %in% colnames(dt))
    stopifnot(id %in% colnames(dt))
    mapping <- dt[[what]]
    names(mapping) <- dt[[id]]
    return(mapping)
}

# map by variant
map_by_variant <- function(variants, mapping, allow_dups=TRUE){
    stopifnot(any(names(mapping) %in% variants))
    mapped <- unlist(lapply(variants, function(v){
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
    
    dt <- fread(input_path)
    annotation <- fread(vep_path)

    # check that all variants are contained
    variants_in_dt <- unlist(strsplit(dt$varid, split = ";"))
    ok <- as.logical(sum(variants_in_dt %in% annotation$varid)/length(variants_in_dt))
    stopifnot(ok)
    
    # create mapping
    map_mac <- create_mapping(annotation, "MAC")
    map_af <- create_mapping(annotation, "info.AF")
    map_maf <- create_mapping(annotation, "MAF")
    map_enstid <- create_mapping(annotation, "csqs.transcript_id")
    map_exon <- create_mapping(annotation, "csqs.exon")
    map_intron <- create_mapping(annotation, "csqs.intron")
    map_revel <- create_mapping(annotation, "csqs.revel_score")
    map_cadd <- create_mapping(annotation, "csqs.cadd_phred")
    map_function_csqs <- create_mapping(annotation, "csqs.most_severe_consequence")
    map_csqs <- create_mapping(annotation, "consequence_category")
    map_codons <- create_mapping(annotation, "csqs.codons")
    map_amino_acids <- create_mapping(annotation, "csqs.amino_acids")

    # map items from vep
    dt$consequence_category <- map_by_variant(dt$varid, map_csqs)
    dt$most_severe_consequence <- map_by_variant(dt$varid, map_function_csqs)
    dt$AF <- map_by_variant(dt$varid, map_af)
    dt$MAF <- map_by_variant(dt$varid, map_maf)
    dt$MAC <- map_by_variant(dt$varid, map_mac)
    dt$transcript <- map_by_variant(dt$varid, map_enstid, allow_dups = FALSE)
    dt$exon <- map_by_variant(dt$varid, map_exon)
    dt$intron <- map_by_variant(dt$varid, map_intron)
    dt$codons <- map_by_variant(dt$varid, map_codons)
    dt$amino_acids <- map_by_variant(dt$varid, map_amino_acids)

    # make columns nicer
    new_order <- c("s","gene_id","transcript","varid","gts","hom_alt_n", "phased.a1", "phased.a2",
                   "unphased.n", "pKO" ,"knockout", "consequence_category", "most_severe_consequence",
                   'AF', 'MAF', 'MAC', 'transcript','exon','intron','codons','amino_acids')
    stopifnot(all(new_order %in% colnames(dt)))
    stopifnot(all(colnames(dt) %in% new_order))
    
    # create subsets and export
    only_plofs <- grepl("pLoF", dt$consequence_category) & !grepl("damaging_missense", dt$consequence_category)
    only_missense <- !grepl("pLoF", dt$consequence_category) & grepl("damaging_missense", dt$consequence_category)
   
    # export all 
    outfile <- paste0(out_prefix,".pLoF_damaging_missense.txt.gz")
    fwrite(dt, outfile)

    # export plofs
    outfile <- paste0(out_prefix,".pLoF.txt.gz")
    fwrite(dt[only_plofs,], outfile)

    # export plofs
    outfile <- paste0(out_prefix,".damaging_missense.txt.gz")
    fwrite(dt[only_missense,], outfile)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "?")
parser$add_argument("--vep_path", default=NULL, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

