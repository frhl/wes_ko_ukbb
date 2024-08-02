# Author: Samvida S. Venkatesh
# Date: 04/05/22

library(argparse)
library(survival)
library(survminer)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
theme_set(theme_bw())

# Parse arguments ----

parser <- ArgumentParser()
parser$add_argument("--geneNames", required=TRUE, 
                    help = "List of genes to test")
parser$add_argument("--chunk", required=TRUE,
                    help = "Results chunk number")
parser$add_argument("--includePRS", required=TRUE,
                    help = "Whether to include PRS as covariate") 
parser$add_argument("--refGroup", required=TRUE,
                    help = "Reference group against which to calculate hazard ratios") 
args <- parser$parse_args()
print(args)

includePRS <- toupper(args$includePRS) %in% c("T", "TRUE")

# Read and wrangle data into correct formats ----

#input_prefix <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/survival/knockouts/pLoF_damaging_missense/"

input_prefix <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/survival/samvida/UKB.carrier_matrix.eur.unrel.af05.pp0.90.pLoF_damaging_missense.keep.wes200k.replication.b1of1."
out_prefix <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/for_fred_ko_project/2311_analyses/results/ref_group_", 
                     args$refGroup, "/unconditioned/chunk", args$chunk)
if (includePRS) {
  out_prefix <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/for_fred_ko_project/2311_analyses/results/ref_group_", 
                       args$refGroup, "/prs_conditioned/chunk", args$chunk)
  prs_prefix <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/prs/scores_new/"
  prs_diags_include <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/prs/validation/ldsc_summary_keep_phenos.txt",
                                  sep = "\t", header = F, stringsAsFactors = F)$V1
}
dir.create(out_prefix)

## General covariates file for Cox regression covariates ----

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

QCOVARS <- strsplit("PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10", 
                    "[+ ]|[+]|\n")[[1]]
CATCOVARS <- strsplit("sex+UKB_assmt_centre+birth_cohort", 
                      "[+ ]|[+]|\n")[[1]]
COX_COVARS <- c(QCOVARS, CATCOVARS)
if (includePRS) COX_COVARS <- c(COX_COVARS, "off_chrom_PRS")
print(COX_COVARS)

## Gene-phenotype combinations to test ----

gene_pheno_tests <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/for_fred_ko_project/2311_analyses/gene_phenos_to_test_ref_group_",
                                      args$refGroup, ".txt"),
                               sep = "\t", header = T, stringsAsFactors = F, quote = "", 
                               comment.char = "$")

gene_pheno_tests <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/for_fred_ko_project/2311_analyses/gene_pheno_counts.txt",
                               sep = "\t", header = T, stringsAsFactors = F, quote = "", 
                               comment.char = "$")

GENES <- strsplit(args$geneNames, "[|]|[| ]|\n")[[1]]
gene_pheno_tests <- subset(gene_pheno_tests, gene_pheno_tests$gene %in% GENES)

DIAGNOSES <- unique(gene_pheno_tests$diagnosis)
gene_pheno_tests <- split(gene_pheno_tests, f = gene_pheno_tests$gene)

## Survival data ----
surv_dat <- 
  readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_survival_dat_all_phenotypes.rds")

# Subset to diagnoses of interest
surv_dat <- surv_dat[DIAGNOSES]

## IDs and classification ----

if (args$refGroup == "het") {
  KOLEVELS <- c("het", "chet_cis", "knockout", "wt")
} else if (args$refGroup == "chet_cis") {
  KOLEVELS <- c("chet_cis", "het", "hom", "chet_trans", "wt")
}

id_class <- lapply(GENES, function (g) {
  # fname <- gzfile(paste0(input_prefix, gene_pheno_tests[[g]]$input_dir[1]), "rt")
  res <- read.table(paste0(input_prefix, g, ".tsv"), 
                    header = T, stringsAsFactors = F,
                    sep = "\t", na.strings = "")[, c("eid", "annotation")]
  res$annotation[res$annotation == "chet"] <- "chet_trans"
  res$annotation[res$annotation == "cis"] <- "chet_cis"
  
  res <- res %>% 
    rename(group = annotation) 
  
  if (args$refGroup == "het") {
    res <- res %>% mutate(group = ifelse(group == "hom" | group == "chet_trans", 
                                         "knockout", group))
  }
  
  res <- res %>% mutate(eid = as.character(eid),
                        group = factor(as.character(group),
                                       levels = KOLEVELS))
  return (res)
})
names(id_class) <- GENES

# Function to add off-chromosome PRS for a given gene-phenotype combination ----

# Read dictionary (to know which phenotypes are Spiros phenos)
spiros_dict <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary_with_unix_codes.txt",
                          sep = "\t", header = T, stringsAsFactors = F, 
                          quote = "", comment.char = "$")

source("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/scripts/post_hoc/utils.R")
mapping <- get_mapping_ensembl_to_contig()

getPRS <- function (gene, prs_fname) {
  prs_fname <- gzfile(prs_fname, "rt")
  prs_dat <- read.table(prs_fname, header = T, stringsAsFactors = F,
                        sep = "\t")
  colnames(prs_dat)[1] <- "eid"
  
  # get chromosome for given gene 
  gchr <- paste0("chr", mapping[gene])
  # remove this column from PRS calculation
  prs_calc <- prs_dat %>% dplyr::select(-all_of(c("eid", gchr)))
  prs_calc <- rowSums(prs_calc)
  
  # Create eid-offchrPRS table
  res <- data.frame(eid = prs_dat$eid,
                    off_chrom_PRS = prs_calc)
  res$eid <- as.character(res$eid)
  res$off_chrom_PRS <- as.numeric(res$off_chrom_PRS)
  return (res)
}

# Add grouping and covariates to survival data and log individuals retained ----

# If there is a requested regression covariate not already in the id files,
# add from general covars (if available)
add_covars <- COX_COVARS[!COX_COVARS %in% colnames(id_class[[1]])]

if ("birth_cohort" %in% add_covars) {
  general_covars <- general_covars %>% 
    mutate(birth_cohort = factor(paste0("<", 
                                        plyr::round_any(year_of_birth, 
                                                        10, f = floor))))
}
add_gen <- general_covars %>% select(any_of(c("eid", add_covars)))

# Ensure quantitative and categorical covars are the right data types
full_id_dat <- lapply(id_class, function (g_df) {
  res <- left_join(g_df, add_gen, by = "eid") %>%
    mutate(across(any_of(QCOVARS), as.numeric)) %>%
    mutate(across(any_of(CATCOVARS), as.factor))
  res <- res[complete.cases(res), ]
  return (res)
})

# Function to run Cox-PH regression ----

COXMODFORM <- paste0("Surv(age_at_first_record, age_at_event, event_observed) ~ group + ",
                     paste0(COX_COVARS, collapse = " + "))

runCoxRegression <- function (gene, diag, prs_fname) {
  df <- surv_dat[[diag]]
  df$eid <- as.character(df$eid)
  df <- inner_join(df, full_id_dat[[gene]], by = "eid")
  
  # Add PRS column if needed
  if (prs_fname != "") {
    df <- inner_join(df, getPRS(gene, prs_fname), by = "eid")
  }
  
  res <- tryCatch({
    cox_mod <- coxph(formula(COXMODFORM), data = df)
    # Get text (HRs and p-values) from Cox regression 
    res_table <- broom::tidy(cox_mod, exp = T) %>%
      mutate(add_text = paste0(gsub("group", "", term), ":",
                               " HR=", signif(estimate, 3), 
                               ", P=", signif(p.value, 3)))
    # Also return table to print
    table_write <- res_table %>%
      filter(grepl("^group", term)) %>%
      mutate(term = gsub("group", "", term),
             gene = gene, diagnosis = diag) %>%
      select(all_of(c("gene", "diagnosis", "term", "estimate", "std.error", "p.value")))
    
    return (table_write)
  }, error = function (cond) {
    return (NULL)
  })
  return (res)
}

# Run all analyses ----

full_res <- lapply(GENES, function (g) {
  # Get list of diagnoses for which to run Cox regression 
  run_diags <- unique(unlist(gene_pheno_tests[[g]]$diagnosis))
  
  table_write <- lapply(run_diags, function (diag) {
    prs_fname <- ""
    if (includePRS) {
      # Only run if the PGS file exists
      # add "spiro_" to diagnosis name if it is not a duncan pheno
      # add combined if not
      if (diag %in% spiros_dict$unix_code) {
        prs_diag <- paste0("spiro_", diag)
      } else {
        prs_diag <- paste0(diag, "_combined")
      }
      # get PRS file for the diagnosis (if it exists)
      if (prs_diag %in% prs_diags_include) {
        prs_fname <- paste0(prs_prefix, prs_diag, "_pgs_chrom.txt.gz")
      }
    }
    coxres <- runCoxRegression(g, diag, prs_fname)
    return (coxres)
  })
  table_write <- bind_rows(table_write)
  return (table_write)
})
full_res <- bind_rows(full_res)
write.table(full_res, paste0(out_prefix, 
                             "/coxph_results.txt"),
            sep = "\t", row.names = F, quote = F)


