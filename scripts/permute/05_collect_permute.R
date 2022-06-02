#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){

    #print(args)  
    stopifnot(dir.exists(args$target_dir))

    # subset to files that are in the main chr folders of target dir
    files <- list.files(args$target_dir, recursive = TRUE, full.names = TRUE, pattern = ".pvalues")
    target_dirs <- list.files(args$target_dir, full.names = TRUE)
    target_regex <- paste(paste0("(",target_dirs[grepl("chr[0-9]+$", target_dirs)],")"), collapse = '|')
    files <- files[grepl(target_regex , files)]

    # iterate over file and aggregate numbers
    final <- do.call(rbind, lapply(files, function(f){
      
      # house keeping 
      main <- stringr::str_extract(basename(f), "ENS.+\\.pvalues\\.txt.gz")
      chrom <- stringr::str_extract(basename(f), "chr[0-9]+")
      gene <- stringr::str_extract(basename(f), "ENSG[0-9]+")
      pheno <- unlist(lapply(strsplit(gsub("\\.pvalues\\.txt\\.gz","",main), split = '\\_'), function(x) paste0(x[-1], collapse = "_")))
      
      # actual empirical P-values
      d <- fread(cmd = paste("zcat", f))
      bool_p_values <- as.logical(d$is_permuted)
      p_true <- d$p[!d$is_permuted]
      p_values <- d$p[bool_p_values]
      p_mean <- mean(p_values)
      p_median <- median(p_values)
      n_permuted <- sum(d$is_permuted)
      n_not_permuted <- sum(!d$is_permuted)
      
      # number of P-values above/below thresholds
      sum_p_values_le_p_true <- sum(p_values <= p_true)
      sum_p_values_lt_p_true <- sum(p_values < p_true)
      sum_p_values_ge_p_true <- sum(p_values >= p_true)
      sum_p_values_gt_p_true <- sum(p_values > p_true)
      sum_p_equal_zero <- sum(p_values == 0)
      sum_p_equal_one <- sum(p_values == 1)
      
      # get deciles
      sum_p_decile_0_1 <- sum(p_values <= 0.0 & p_values < 0.1)
      sum_p_decile_1_2 <- sum(p_values <= 0.1 & p_values < 0.2)
      sum_p_decile_2_3 <- sum(p_values <= 0.2 & p_values < 0.3)
      sum_p_decile_3_4 <- sum(p_values <= 0.3 & p_values < 0.4)
      sum_p_decile_4_5 <- sum(p_values <= 0.4 & p_values < 0.5)
      sum_p_decile_5_6 <- sum(p_values <= 0.5 & p_values < 0.6)
      sum_p_decile_7_8 <- sum(p_values <= 0.7 & p_values < 0.8)
      sum_p_decile_8_9 <- sum(p_values <= 0.8 & p_values < 0.9)
      sum_p_decile_9_10 <- sum(p_values < 0.9 & p_values <= 1.0)
    
      # summarize in out data
      d_out <- data.frame(
          chrom,
          gene,
          pheno,
          p_true,
          p_mean,
          p_median,
          n_permuted,
          n_not_permuted,
          sum_p_values_le_p_true,
          sum_p_values_lt_p_true,
          sum_p_values_ge_p_true,
          sum_p_values_gt_p_true,
          sum_p_equal_zero,
          sum_p_equal_one,
          sum_p_decile_0_1,
          sum_p_decile_1_2,
          sum_p_decile_2_3,
          sum_p_decile_3_4,
          sum_p_decile_4_5,
          sum_p_decile_5_6,
          sum_p_decile_7_8,
          sum_p_decile_8_9,
          sum_p_decile_9_10,
          files = f
      )
      
      return(d_out)
    })) 
  
  outfile <- paste0(args$out_prefix, ".txt.gz")
  fwrite(final, outfile, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--target_dir", default=NULL, help = "chromosome")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

