# Author: Samvida S. Venkatesh
# Date: 04/05/2022

library(tidyverse)

#ALL_GENES <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/for_fred_ko_project/2311_analyses/genes_to_test.txt",
ALL_GENES <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/survival/2401_cox_replication_genes.txt",
                        header = F, stringsAsFactors = F, sep = "\t")$V1

#submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/for_fred_ko_project/scripts/submit_coxph_models.sh"
submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/scripts/survival/replicate_missing/submit_coxph_models.sh"

# Types of files
#ftypes <- c("remove.wes200k.replication_full",
#            "remove.wes200k.replication",
#            "keep.wes200k.replication_full",
#            "keep.wes200k.replication")

ftypes <- c("remove.wes200k",
            "keep.wes200k")


# Reference groups
refgps <- c("het", "chet_cis")

# Split into chunks of N genes
NGENE_CHUNK <- 16
gene_chunks <- split(ALL_GENES, ceiling(seq_along(ALL_GENES)/NGENE_CHUNK))

for (ref_gp in refgps) {
  for (filetype in ftypes) {
    #for (i in c(1)){
    for (i in 1:length(gene_chunks)) {
      genes_in_chunk <- paste0(gene_chunks[[i]], collapse = "|")
      # Submit job parameters
      job_options <- paste0(
        "--export=",
        paste0(
          "geneNames=\"", genes_in_chunk, "\",",
          "chunk=\"", i, "\",",
          "fileType=\"", filetype, "\",",
          "refGroup=\"", ref_gp, "\""
        )
      )
      job_submission <- paste("sbatch", job_options, submission_script)
      system(job_submission)
      print(job_submission)
    }
  }
}



