
# 
devtools::load_all("utils/modules/R/gwastools")
library(data.table)
library(argparse)

# map from ENSEMBL to HGNC
bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensembl_to_hgnc <- bridge$hgnc_symbol
names(ensembl_to_hgnc) <- bridge$ensembl_gene_id
ensembl_to_pos <- (bridge$start_position + bridge$end_position)/2
names(ensembl_to_pos) <- bridge$ensembl_gene_id
ensembl_to_contig <- bridge$chromosome_name
names(ensembl_to_contig) <- bridge$ensembl_gene_id

process_saige_df <- function(d, f){
    # deal with conditional P-values
    if ("p.value_c" %in% colnames(d)) {
        d$p.value <- d$p.value_c
        d$p.value.NA <- d$p.value.NA_c
        d$BETA <- d$BETA_c
        d$SE <- d$SE_c
        d$Tstat <- d$Tstat_c
        d$var <- d$var_c
    }
    
    # 
    d$p.value <- as.numeric(d$p.value)
    d$p.value.NA <- as.numeric(d$p.value.NA)
    d$BETA <- as.numeric(d$BETA)
    d$SE <- as.numeric(d$SE)

    d$ensembl_gene_id <- d$MarkerID
    d$ensembl_gene_id[!grepl("ENSG", d$ensembl_gene_id)] <- NA
    d$pvalue = d$p.value
    d$prs <- ifelse(grepl("locoprs", f),"With PRS", "Without PRS")
    d$analysis <- basename(f)
    d$analysis <- stringr::str_extract(d$analysis, "200k_.+pLoF_damaging_missense")
    d$analysis <- gsub("200k_", "", d$analysis)
    d$analysis <- gsub("_pLoF_damaging_missense", "", d$analysis)
    d$phenotype <- d$analysis

    # add hgnc and grch38 position
    d$hgnc_symbol <- ensembl_to_hgnc[d$ensembl_gene_id]
    d$contig <- ensembl_to_contig[d$ensembl_gene_id]
    d$pos <- ensembl_to_pos[d$ensembl_gene_id]
    d$N_ko_case <- d$N_case_hom #unlist(ifelse("N_ko_case" %in% colnames(d), list(d$N_case_hom), NA))
    d$N_ko_ctrl <- d$N_ctrl_hom #unlist(ifelse("N_ko_ctrl" %in% colnames(d), list(d$N_ctrl_hom), NA))
    d$N_ko <- d$AC_Allele2/2
    d$KF_ko <- d$AF_Allele2
    d$KF_case <- d$AF_case
    d$KF_ctrl <- d$AF_ctrl
    return(d)
}


get_formatted_df <- function(files){
    d <- do.call(rbind, lapply(files, function(f){
        stopifnot(file.exists(f))
        d <- fread(f)
        if (is.numeric(d$p.value)){
            if (nrow(d) > 0){
                return(process_saige_df(d, f))
            } else {
                return(NULL)
            }
        } else {
            write(paste0(f, "was excluded because p-value col was invalid."), stderr())   
        }
    }))
    return(d)
}


main <- function(args){

    N_ko_cutoff <- as.numeric(args$N_ko_cutoff)
    N_ko_case_cutoff <- as.numeric(args$N_ko_case_cutoff)
    p_cutoff <- args$p_cutoff

    # ensure we are using saige with PRS when possible
    ref <- fread(args$ref_file)
    ref$ID <- paste0(ref$MarkerID,"_",ref$phenotype)
    mapping_n_ko_case <- ref$N_ko_case
    mapping_n_ko_ctrl <- ref$N_ko_ctrl
    names(mapping_n_ko_case) <- ref$ID
    names(mapping_n_ko_ctrl) <- ref$ID 

    # we only include things that were significant in the primary analysis
    d <- fread(args$merged_hits)
    d$ID <- paste0(d$MarkerID,"_",d$phenotype)
    d$N_case_hom <- mapping_n_ko_case[d$ID]
    d$N_ctrl_hom <- mapping_n_ko_ctrl[d$ID]
    d <- d[d$ID %in% ref$ID,]     

    # wei has not added N_ko_case and N_ko_ctrl so we need to grab them from the past file
    d <- rbindlist(lapply(1:nrow(d), function(idx){
        dt <- d[idx,]
        path <- dt$filepath
        return(process_saige_df(dt, path))    
    }))

    # perform subsets
    d <- d[d$N_ko >= N_ko_cutoff,]
    d <- d[d$N_ko_case >= N_ko_case_cutoff,]

    # p-value cutoff
    if (!is.null(p_cutoff)){
        p_cutoff <- as.numeric(p_cutoff)
        d <- d[(d$p.value < p_cutoff),]
    }

    # columsn to keep
    cols <- c(
        'phenotype','CHR','MarkerID','hgnc_symbol',
        'MissingRate','BETA','SE','Tstat',
        'var','p.value','p.value.NA','N_case',
        'N_ctrl','N_ko_case','N_ko_ctrl','N_ko',
        'KF_ko','KF_case','KF_ctrl','prs')
    
    # subset and order
    d <- d[,..cols]
    d <- d[order(d$p.value),]
    d$hgnc_symbol[d$hgnc_symbol==""] <- NA
    d$hgnc_symbol[d$hgnc_symbol==" "] <- NA
    outfile <- paste0(args$out_prefix, ".txt.gz")
    write(paste0("writing ", outfile), stdout())
    fwrite(d, outfile, sep = "\t", na = "n/a")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--ref_file", default=NULL, required = TRUE, help = "")
parser$add_argument("--merged_hits", default=NULL, required = TRUE, help = "")
parser$add_argument("--N_ko_case_cutoff", default=NULL, required = TRUE, help = "")
parser$add_argument("--N_ko_cutoff", default=NULL, required = TRUE, help = "")
parser$add_argument("--p_cutoff", default=NULL, required = FALSE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)


