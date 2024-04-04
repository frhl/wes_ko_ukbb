
library(tidyverse)
source("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/scripts/post_hoc/utils.R")

# read spiro unix_encoding
encoding <- fread("data/phenotypes/phenotype_icd_chapter.txt")

# Read results tables ----
uncond <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/for_fred_ko_project/2301_analyses/results/ref_group_het/unconditioned/all_coxph_results.txt",
                  sep = "\t", header = T, stringsAsFactors = F, quote = "", 
                  comment.char = "$")
uncond <- uncond %>% 
  mutate(test_combo = paste0(gene, ":", diagnosis),
         analysis_type = "unconditioned")

cond <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/for_fred_ko_project/2301_analyses/results/ref_group_het/prs_conditioned/all_coxph_results.txt",
                     sep = "\t", header = T, stringsAsFactors = F, quote = "", 
                     comment.char = "$")
cond <- cond %>% 
  mutate(test_combo = paste0(gene, ":", diagnosis),
         analysis_type = "conditioned")

# Combine results
dat <- bind_rows(uncond, cond)
dat <- dat %>% 
  pivot_wider(values_from = c(estimate, std.error, p.value),
              names_from = analysis_type)

mapping <- get_mapping_ensembl_to_hgnc()
dat$hgnc_symbol <- mapping[dat$gene]

# Set multiple-testing thresholds 
unique_comb_tests <- unique(uncond$test_combo)
PTHRESH <- 0.05/length(unique_comb_tests)

# Forest plot for hazard ratios ----

# need to include this to not remove hits that don't have PRS
dat$p.value_conditioned[is.na(dat$p.value_conditioned)] <- dat$p.value_unconditioned[is.na(dat$p.value_conditioned)]
dat$min_p <- pmin(dat$p.value_unconditioned, dat$p.value_conditioned)
combos_plot <- dat$test_combo[which(dat$term == "knockout" & dat$min_p <= PTHRESH)]

# get things to save
to_save <- dat[dat$test_combo %in% combos_plot,]
to_save <- to_save[to_save$term %in% "knockout",]
to_save$MarkerID <- to_save$gene
to_save$phenotype <- to_save$diagnosis
to_save$CHR <- paste0("chr",ensembl_to_contig[to_save$gene])
to_save <- to_save[,c("MarkerID","hgnc_symbol","phenotype", "CHR")]

# perform little mapping
encoding$unix_code_no_spiro <- gsub("spiro_","",encoding$unix_code)
mapping_no_spiro_to_spiro <- encoding$unix_code
names(mapping_no_spiro_to_spiro) <- encoding$unix_code_no_spiro
to_save$spiro_phenotype <- mapping_no_spiro_to_spiro[to_save$phenotype]
bool_na <- is.na(to_save$spiro_phenotype)
to_save$spiro_phenotype[bool_na] <- paste0(to_save$phenotype[bool_na],"_combined")
# remove spiro phenotypes
to_save$phenotype <- to_save$spiro_phenotype
to_save$spiro_phenotype <- NULL

fwrite(to_save, file = "data/survival/results/2302_sig_assoc.txt.gz")

