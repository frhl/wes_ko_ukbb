
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

# icd codes for phenotypes
icd <- fread("data/phenotypes/phenotype_icd_chapter.txt")
icd <- icd[order(icd$ICD_chapter),]
na_chapter <- "None of the above"
icd[is.na(icd$ICD_chapter_desc_short)]$ICD_chapter_desc_short <- na_chapter
chapters <- unique(icd$ICD_chapter_desc_short)

# columsn to keep
cols_rec_encoding <- c(
    'phenotype','CHR','MarkerID','hgnc_symbol',
    'MissingRate','BETA','SE','Tstat',
    'var','p.value','p.value.NA','N_case',
    'N_ctrl','N_ko_case','N_ko_ctrl','N_ko',
    'KF_ko','KF_case','KF_ctrl','prs')

cols_add_encoding <- c(
    'phenotype','CHR','MarkerID','hgnc_symbol',
    'MissingRate','BETA','SE','Tstat',
    'var','p.value','p.value.NA','N_case',
    'AF_Allele2','AC_Allele2',
    'AF_case', 'AF_ctrl', 'N_ctrl', 
    'N_case_hom', 'N_case_het', 'N_ctrl_hom', 'N_ctrl_het',
    'prs')


process_saige_df <- function(d, f, gt_encoding = "001"){
    # deal with conditional P-values
    if ("p.value_c" %in% colnames(d)) {
        write(paste("Detected conditional P.value (p.value_c) for ",f), stderr())
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
    # knockout encoding
    if (gt_encoding == "001"){
        d$N_ko_case <- d$N_case_hom #unlist(ifelse("N_ko_case" %in% colnames(d), list(d$N_case_hom), NA))
        d$N_ko_ctrl <- d$N_ctrl_hom #unlist(ifelse("N_ko_ctrl" %in% colnames(d), list(d$N_ctrl_hom), NA))
        d$N_ko <- d$AC_Allele2/2
        d$KF_ko <- d$AF_Allele2
        d$KF_case <- d$AF_case
        d$KF_ctrl <- d$AF_ctrl
        d <- d[,..cols_rec_encoding]
    } else if (gt_encoding == "012"){
        d <- d[,..cols_add_encoding]
    }
    return(d)
}


get_formatted_df <- function(files, gt_encoding){
    d <- do.call(rbind, lapply(files, function(f){
        stopifnot(file.exists(f))
        d <- fread(f)
        if (is.numeric(d$p.value)){
            if (nrow(d) > 0){
                return(process_saige_df(d, f, gt_encoding))
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
    p_cutoff <- as.numeric(args$p_cutoff)
    path_header <- args$path_header
    gt_encoding <- args$gt_encoding
    header <- fread(path_header, header=FALSE)$V1

    # get files
    files <- gwastools::list_files_saige(cond = args$cond, prs = args$prs)
    files_found <- length(files)
    if (files_found < 2) stop(paste("Only",files_found,"files were found!"))

    # read files and format
    d <- get_formatted_df(files, gt_encoding)
    d <- d[d$phenotype %in% header,]
 
    # perform subsets
    if (gt_encoding == "001") { 
        d <- d[d$N_ko >= N_ko_cutoff,]
        d <- d[d$N_ko_case >= N_ko_case_cutoff,]
    } else if (gt_encoding == "012") {
        d <- d[d$AC_Allele2 >= (N_ko_cutoff*2),]
        d <- d[d$N_case_hom >= N_ko_case_cutoff,]
    }

    # add ICD code
    if (FALSE) {
        d <- merge(icd, d, by.x = "unix_code", by.y = "phenotype", all.y=TRUE)
        d$ICD_chapter_desc[is.na(d$ICD_chapter_desc)] <- na_chapter
        d$ICD_chapter_desc_short[is.na(d$ICD_chapter_desc_short)] <- na_chapter
    }

    # re-order and subset
    d <- d[d$p.value < p_cutoff,]
    d <- d[order(d$p.value), ]

    # clean up
    d$hgnc_symbol[d$hgnc_symbol==""] <- NA
    d$hgnc_symbol[d$hgnc_symbol==" "] <- NA
    outfile <- paste0(args$out_prefix, ".txt.gz")
    write(paste0("writing ", outfile), stdout())
    fwrite(d, outfile, sep = "\t", na = "n/a")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_header", default=NULL, required = TRUE, help = "")
parser$add_argument("--p_cutoff", default=NULL, required = TRUE, help = "")
parser$add_argument("--N_ko_case_cutoff", default=NULL, required = TRUE, help = "")
parser$add_argument("--N_ko_cutoff", default=NULL, required = TRUE, help = "")
parser$add_argument("--cond", default="none", required = TRUE, help = "")
parser$add_argument("--prs", default=NULL, required = TRUE, help = "")
parser$add_argument("--gt_encoding", default="001", required = TRUE, help = "A character string to indicate genotype encoding, either recessive '001' or additive '012'.")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)


