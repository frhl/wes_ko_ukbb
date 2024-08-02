# Author: Samvida S. Venkatesh
# Date: 04/05/22

library(argparse)
library(survival)
#library(survminer)
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
parser$add_argument("--fileType", required=TRUE,
                    help = "File prefix for annotation and sample type") 
parser$add_argument("--refGroup", required=TRUE,
                    help = "Reference group against which to calculate hazard ratios") 
args <- parser$parse_args()
print(args)

# Read and wrangle data into correct formats ----

input_prefix <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/survival/samvida/240226/UKB.carrier_matrix.qced.eur.unrel.af05.pp0.90.pLoF_damaging_missense.", args$fileType, 
                       ".b1of1.")
out_prefix <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/survival/missing_replication/240226_analyses/results/ref_group_", 
                     args$refGroup, "/unconditioned/", args$fileType, "_chunk", args$chunk)

dir.create(dirname(out_prefix), recursive = TRUE)

## General covariates file for Cox regression covariates ----

general_covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220504_QCd_demographic_covariates.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
general_covars$eid <- as.character(general_covars$eid)

QCOVARS <- strsplit("PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10", 
                    "[+ ]|[+]|\n")[[1]]
CATCOVARS <- strsplit("sex+UKB_assmt_centre+birth_cohort", 
                      "[+ ]|[+]|\n")[[1]]
COX_COVARS <- c(QCOVARS, CATCOVARS)
print(COX_COVARS)

## Gene-phenotype combinations to test ----

GENES <- strsplit(args$geneNames, "[|]|[| ]|\n")[[1]]
DIAGNOSES <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/survival/2401_cox_replication_phenos.txt",
                        header = F, stringsAsFactors = F, sep = "\t")$V1

write("Skipping ENSG00000127423..", stderr())
GENES <- GENES[! GENES %in% "ENSG00000127423"]

# standardize by removix refix or suffix
DIAGNOSES  <- gsub("^spiro_", "", gsub("_combined$", "", DIAGNOSES))

## Survival data ----
surv_dat <- 
  readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_survival_dat_all_phenotypes.rds")

# check id diagnoses is not present
stopifnot(DIAGNOSES %in% names(surv_dat)) # error out
surv_dat <- surv_dat[DIAGNOSES]


## IDs and classification ----

if (args$refGroup == "het") {
  KOLEVELS <- c("het", "chet_cis", "knockout", "wt")
} else if (args$refGroup == "chet_cis") {
  KOLEVELS <- c("chet_cis", "het", "hom", "chet_trans", "wt")
}

id_class <- lapply(GENES, function (g) {
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

runCoxRegression <- function (gene, diag) {
  df <- surv_dat[[diag]]
  df$eid <- as.character(df$eid)
  print(head(df))
  print(head(full_id_dat[[gene]]))
  
  df <- inner_join(df, full_id_dat[[gene]], by = "eid")
  
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
  table_write <- lapply(DIAGNOSES, function (diag) {
    coxres <- runCoxRegression(g, diag)
    return (coxres)
  })
  table_write <- bind_rows(table_write)
  return (table_write)
})
full_res <- bind_rows(full_res)
write.table(full_res, paste0(out_prefix, 
                             "_coxph_results.txt"),
            sep = "\t", row.names = F, quote = F)


