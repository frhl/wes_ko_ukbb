
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


list_files_saige <- function(cond="none", prs="include", regex = "\\.txt\\.gz"){

    # set up paths
    wd <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb"
    step2_dir <- file.path(wd, "data/saige/output/binary/step2/min_mac4")
    step2_common_dir <- file.path(wd, "data/saige/output/binary/step2_common/min_mac4")
    step2_rare_dir <- file.path(wd, "data/saige/output/binary/step2_rare_cond/min_mac4")
    step2_combined_dir <- file.path(wd, "data/saige/output/binary/step2_rare_cond/min_mac4")
    #step2_collapsed_dir <- file.path(wd, "data/saige/output/binary/step2_collapsed/min_mac4")
    
    # subset paths based on condions
    if (cond %in% "none"){
        files <- list.files(step2_dir, full.names = TRUE, pattern = regex)
    } else if (cond %in% "common"){
        files <- list.files(step2_common_dir, full.names = TRUE, pattern = regex)
    } else if (cond %in% "rare"){
        files <- list.files(step2_rare_dir, full.names = TRUE, pattern = regex)
    } else if (cond %in% "combined"){
        files <- list.files(step2_combined_dir, full.names = TRUE, pattern = regex)
    } else {
        stop(paste(cond, "is not a valid. Try 'none','common','rare' or 'combined'"))
    }
    
    # perform final subset by PRS
    files <- sort(files)
    is_prs <- grepl("locoprs.txt.gz", files)
    if (prs %in% "include"){
        return(files)
    } else if (prs %in% "exclude"){
        # exclude any PRS files
        return(files[!is_prs])
    } else if (prs %in% "only") {
        # only include PRS
        return(files[is_prs])
    } else if (prs %in% "prefer") {
        return(gwastools::unique_but_prefer_regex(files, regex="_locoprs"))
    } else {
        stop(paste(prs, "is not valid. Must be either 'include','exclude','prefer' or 'only'."))
    }
    
}

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
    p_cutoff <- as.numeric(args$p_cutoff)
    path_header <- args$path_header
    header <- fread(path_header, header=FALSE)$V1

    # get files
    files <- list_files_saige(cond = args$cond, prs = args$prs)
    print(head(files))
    
    # read files and format
    d <- get_formatted_df(files)
    d <- d[d$phenotype %in% header,]
  
    # perform subsets
    d <- d[d$N_ko >= N_ko_cutoff,]
    d <- d[d$N_ko_case >= N_ko_case_cutoff,]
    d <- d[d$p.value < p_cutoff,]
    d <- d[order(d$p.value), ]

    # columsn to keep
    cols <- c(
        'phenotype','CHR','MarkerID','hgnc_symbol',
        'MissingRate','BETA','SE','Tstat',
        'var','p.value','p.value.NA','N_case',
        'N_ctrl','N_ko_case','N_ko_ctrl','N_ko',
        'KF_ko','KF_case','KF_ctrl','prs')
    
    # subset and order
    d <- d[,..cols]

    # get tier A in the top
    d <- d[order(d$p.value),]
    d1 <- d[d$N_ko_case == 1,]
    d2 <- d[d$N_ko_case >= 2,]
    d <- rbind(d2, d1)
        
    # put
    outfile <- paste0(args$out_prefix, ".txt.gz")
    write(paste0("writing ", outfile), stdout())
    fwrite(d, outfile, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_header", default=NULL, required = TRUE, help = "")
parser$add_argument("--p_cutoff", default=NULL, required = TRUE, help = "")
parser$add_argument("--N_ko_case_cutoff", default=NULL, required = TRUE, help = "")
parser$add_argument("--N_ko_cutoff", default=NULL, required = TRUE, help = "")
parser$add_argument("--cond", default="none", required = TRUE, help = "")
parser$add_argument("--prs", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)


