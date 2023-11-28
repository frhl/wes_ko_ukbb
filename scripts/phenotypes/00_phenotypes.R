#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# for performing residualizing
devtools::load_all('utils/modules/R/gwastools')


# load spiro phenotypes based on Samvida's paths
get_spiro_phenos <- function(){
    
    write("Loading spiros phenotypes..", stderr())
      
    # setup paths
    path_phenotype_matrix <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_phenotype_matrix.txt"
    path_disease_dict <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt"
    stopifnot(file.exists(path_phenotype_matrix))
    stopifnot(file.exists(path_disease_dict))
    
    # load main phenotype matrix
    eid_pheno_matrix <- read.table(
        path_phenotype_matrix,
        sep = "\t", 
        header = T,
        stringsAsFactors = F
    )
    colnames(eid_pheno_matrix) <- gsub("^X", "", colnames(eid_pheno_matrix))
    
    # Disease dictionary
    dictionary <- read.table(
        path_disease_dict,
        sep = "\t", 
        header = T, 
        stringsAsFactors = F,
        quote = "", 
        comment.char = "$")
    PHENOTYPES <- dictionary$phenotype[match(colnames(eid_pheno_matrix)[-1], dictionary$unique_code)]
    colnames(eid_pheno_matrix)[-1] <- PHENOTYPES
    return(eid_pheno_matrix)
    
}

# load brava phenotypes based on Duncan's paths
get_brava_phenos <- function(){
    
    write("Loading brava phenotypes..", stderr())
    
    # paths
    path_phenotypes <- "/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/BRaVa_phenotypes_with_superpopulation_labels.tsv"
    eid_pheno_matrix <- fread(path_phenotypes)
    
    # remove columns not needed
    eid_pheno_matrix$ICD9_string <- NULL
    eid_pheno_matrix$ICD10_string <- NULL
    
    # we remove everything after classifcation column
    idx_keep <- 1:which(grepl('classification', colnames(eid_pheno_matrix)))[1]-1
    eid_pheno_matrix <- eid_pheno_matrix[,..idx_keep]
   
    # this is a duplicate after cleaning headers
    target <- "Female infertility (anatomic causes only)" 
    new <- "female_infertility_anatomic_causes_only"
    colnames(eid_pheno_matrix)[colnames(eid_pheno_matrix) == target] <- new
    
    return(eid_pheno_matrix)
    
}


# clean header
to_unix_friendly_colnames <- function(names){
    cols <- gsub("\\'", "", names)
    cols <- gsub("\\(.+\\)","", cols) 
    cols <- gsub('\\"', "", cols)
    cols <- gsub("\\-","_", cols)
    cols <- gsub("\\&", "_and_", cols)
    cols <- gsub("\\,","_",cols)
    cols <- gsub("(^\\ *)|(\\ *$)","", cols)
    cols <- gsub("\\/", "_", cols)
    # preserve unicode charcaters
    cols <- stringi::stri_trans_general(cols, "latin-ascii")
    # clean up spaces
    cols <- gsub("\\ +", " ", cols)
    cols <- gsub("\\ ", "_", cols)
    cols <- gsub("\\_+", "_", cols)
    cols <- tolower(cols)
    return(cols)
}


main <- function(args){
 
    print(args)  
    stopifnot(file.exists(args$input_path))
    stopifnot(dir.exists(dirname(args$out_path)))
    stopifnot(!is.null(args$covariates))

    # read in data 
    dt <- fread(args$input_path)
    
    # Include spiro phenotypes
    if (args$include_spiros) {
        spiro <- as.data.frame(get_spiro_phenos())
        spiro[,-1] <- spiro[,-1] > 0
        colnames(spiro)[-1] <- paste0("spiro_", to_unix_friendly_colnames(colnames(spiro)))[-1]
        dt <- merge(dt, spiro, by = "eid", all.x = TRUE)
        print(head(dt))
    }

    # include BRAVA pgenotypes
    if (args$include_brava) {
        brava <- as.data.frame(get_brava_phenos())
        brava[,-1] <- brava[,-1] > 0
        colnames(brava)[-1] <- paste0("brava_", to_unix_friendly_colnames(colnames(brava)))[-1]
        dt <- merge(dt, brava, by = "eid", all.x = TRUE)
        print(head(dt))
    }

    # include BRAVA pgenotypes
    if (args$include_lindgren) {
        stopifnot(file.exists(args$input_path_cts_to_bin))
        lindgren <- fread(args$input_path_cts_to_bin)
        phenos <- c("BMI","WHR")
        cols_keep <- colnames(lindgren) %in% c("eid","sex", phenos)
        lindgren_phenos <- lindgren[,cols_keep,with=FALSE]
        # define BMI cases
        lindgren_phenos$lindgren_obesity <- NA
        lindgren_phenos$lindgren_obesity[lindgren_phenos$BMI >= 35] <- TRUE
        lindgren_phenos$lindgren_obesity[(lindgren_phenos$BMI >= 18.5) & (lindgren_phenos$BMI < 25)] <- FALSE
        # define WHR cases
        lindgren_phenos$lindgren_whr <- NA
        lindgren_phenos$lindgren_whr[(lindgren_phenos$WHR >= 0.90) & (lindgren_phenos$sex == 1) ] <- TRUE # men case
        lindgren_phenos$lindgren_whr[(lindgren_phenos$WHR >= 0.85) & (lindgren_phenos$sex == 0) ] <- TRUE # woman case
        lindgren_phenos$lindgren_whr[(lindgren_phenos$WHR < 0.90) & (lindgren_phenos$sex == 1) ] <- FALSE # man ctrl
        lindgren_phenos$lindgren_whr[(lindgren_phenos$WHR < 0.80) & (lindgren_phenos$sex == 0) ] <- FALSE # woman ctrl
        lindgren_phenos <- lindgren_phenos[,c(1,5,6)]
        dt <- merge(dt, lindgren_phenos, by = "eid", all.x = TRUE)
        print(head(dt)) 
    }

    if (!is.null(args$case_count_cutoff)) {    
        
        write("Removing phenotypes by case count cutoff", stderr())  

        # drop phenotypes based on case/control count
        pheno_cols <- names(which(sapply(dt, is.logical))) 
        n_cutoff <- as.numeric(args$case_count_cutoff)

        # load samples that are QCed
        qc_samples <- fread(args$qc_samples)$V1
        

        # discard phenotypes with less than 50 cases in quality controlled samples (200k)
        dt_subset_by_qced_samples <- dt[dt$eid %in% qc_samples,]
        lst <- lapply(pheno_cols, function(ph) sum(dt_subset_by_qced_samples[[ph]], na.rm = TRUE))
        names(lst) <- pheno_cols
        dlst <- stack(lst)
        colnames(dlst) <- c("cases", "phenotype")
        bool_phenos_discard <- dlst$cases < n_cutoff
        dlst <- dlst[!bool_phenos_discard, ]
        write(paste("Discarded", sum(bool_phenos_discard), "phenotypes with less than",n_cutoff,"cases."), stdout())
        phenos_keep <- dlst$phenotype 
        
        # final filtering         
        cols_keep <- (!(colnames(dt) %in% pheno_cols)) | (colnames(dt) %in% phenos_keep)
        dt <- dt[,..cols_keep]
    }
    
    # append new covariates
    dt$age2 <- dt$age ^ 2
    dt$age3 <- dt$age ^ 3
    dt$sex_age <- dt$age * ifelse(dt$sex == 1, 1, -1)

    # only keep samples in which we have no missing covariates 
    covars <- unlist(strsplit(args$covariates, split = ","))
    missing <- rowSums(do.call(cbind, lapply(covars, function(col) is.na(dt[[col]]) ))) > 0
    n_missing <- sum(missing)
    if ((n_missing > 0) & (args$row_na_action == "remove")) {
        write(paste("Note: dropping", n_missing, "sample(s) with missing covariates."),stderr())
        dt <- dt[!missing, ]
    }

    # transform phenotypes (either RINT or INT)
    if (!is.null(args$transform_method) & !is.null(args$transform)){
       
        # only transform selected phenotypes 
        phenotypes <- unlist(strsplit(args$transform, split = ","))
        phenotypes <- phenotypes[phenotypes %in% colnames(dt)]
        stopifnot(length(phenotypes) > 0)

        # do transformation and aggregate
        f <- ifelse(args$transform_method == "int", gwastools::get_int, gwastools::get_rint)
        transformed_phenos <- lapply(phenotypes, function(ph){ return(f(dt[[ph]])) })
        transformed_phenos <- do.call(cbind, transformed_phenos)
        colnames(transformed_phenos) <- paste0(phenotypes, "_",args$transform_method)
        dt <- cbind(dt, transformed_phenos)
    }

    # Finally write file
    nrow_dt <- nrow(dt)
    ncol_dt <- ncol(dt)
    size_msg <- paste(nrow_dt, "samples x ",ncol_dt,"phenotypes")
    write(paste0("Done! Writing ", size_msg, " to ",args$out_path), stdout())
    fwrite(dt, args$out_path, sep = "\t", quote = FALSE)
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_path", default=NULL, help = "Where should the file be written")
parser$add_argument("--input_path", default=NULL, help = "Curated phenotypes input path")
parser$add_argument("--input_path_cts_to_bin", default=NULL, help = "Path to cts phenotypes that should be made into binary (See notion 22-12-14)")
parser$add_argument("--transform_method", default=NULL, help = "Transformation of residuals, either INT or RINT")
parser$add_argument("--transform", default=NULL, help = "What phenotypes should be transformed?")
parser$add_argument("--row_na_action", default="keep", help = "set to 'remove' to delete rows with missing covariates")
parser$add_argument("--covariates", default=NULL, help = "comma seperated string of covariates")
parser$add_argument("--qc_samples", default=NULL, help = "link to duncan's quality control samples (only used to subset by phenotype case count!)")
parser$add_argument("--case_count_cutoff", default=NULL, help = "Discard binary phenotypes with less than n cases.")
parser$add_argument("--include_spiros", default=FALSE, action='store_true', help = "Include spiros phenotypes.")
parser$add_argument("--include_brava", default=FALSE, action='store_true', help = "Include Brava phenotypes")
parser$add_argument("--include_lindgren", default=FALSE, action='store_true', help = "Include cecilia's phenotypes")
args <- parser$parse_args()

main(args)

