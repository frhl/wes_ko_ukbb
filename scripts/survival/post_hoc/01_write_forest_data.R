library(data.table)
library(argparse)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
theme_set(theme_bw())

# get mapping to HGNC
source("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/scripts/post_hoc/utils.R")
mapping <- get_mapping_ensembl_to_hgnc()


# adjust p-values by n number of tests.
p.adjust.n <- function(x, n=NULL, ...){
    if (is.null(n)) return(p.adjust(x, ...))
    if (n < length(x)) stop("param 'n' must be larger than length of x")
    n_to_append <- as.integer(n-length(x))
    to_append <- rep(1, n_to_append)
    new_x <- c(x, to_append)
    q <- p.adjust(new_x, ...)
    q <- q[1:length(x)]
    return(q)
}


main <- function(args){

    in_dir_cond <- args$in_dir_cond
    in_dir_uncond <- args$in_dir_uncond
    out_prefix <- args$out_prefix
    n_tests <- as.numeric(args$n_tests)
    export <- args$export
    stopifnot(n_tests>0)

    stopifnot(dir.exists(in_dir_cond))
    stopifnot(dir.exists(in_dir_uncond))
    
    cond_files <- list.files(recursive = TRUE, full.names=TRUE, in_dir_cond)
    uncond_files <- list.files(recursive = TRUE, full.names=TRUE, in_dir_uncond)

    print(in_dir_cond)
    print(in_dir_uncond)

    print(head(cond_files))
    print(head(uncond_files))


    cond <- rbindlist(lapply(cond_files, function(f) read.table(f, sep = "\t", header = T, stringsAsFactors = F, quote = "", comment.char = "$")))
    uncond <- rbindlist(lapply(uncond_files, function(f) read.table(f, sep = "\t", header = T, stringsAsFactors = F, quote = "", comment.char = "$")))

    # extract files based on PRS conditioning
    uncond <- uncond %>% 
        mutate(test_combo = paste0(gene, ":", diagnosis),
             analysis_type = "unconditioned")
                             
    cond <- cond %>% 
      mutate(test_combo = paste0(gene, ":", diagnosis),
             analysis_type = "conditioned")

    uncond$FDR <- p.adjust.n(uncond$p.value, n=n_tests, method="fdr")
    uncond$FDR[uncond$term != "knockout"] <- NA
    cond$FDR <- p.adjust.n(cond$p.value, n=n_tests, method="fdr")
    cond$FDR[cond$term != "knockout"] <- NA

    # Combine results
    dat <- bind_rows(uncond, cond)
    dat <- dat %>% 
      pivot_wider(values_from = c(estimate, std.error, p.value, FDR),
                  names_from = analysis_type)

    print(head(as.data.frame(dat[dat$diagnosis == "cataract",])))

    # Set multiple-testing thresholds 
    PTHRESH <-  0.05 / n_tests  

    # Forest plot for hazard ratios ----
    dat$hgnc_symbol <- mapping[dat$gene]

    # need to include this to not remove hits that don't have PRS
    dat$p.value_conditioned[is.na(dat$p.value_conditioned)] <- dat$p.value_unconditioned[is.na(dat$p.value_conditioned)]
    dat$min_p <- pmin(dat$p.value_unconditioned, dat$p.value_conditioned)
    dat$FDR <- p.adjust.n(dat$p.value_conditioned, n=n_tests, method="fdr")
    dat$FDR[dat$term != "knockout"] <- NA

    # export significant things or everything
    if (is.null(export)) {
        combos_plot <- dat$test_combo[
            which(
                (dat$term == "knockout") & 
                ((dat$min_p <= PTHRESH) | (dat$FDR < 0.05))
            )]
    } else if (export %in% "everything") {
        combos_plot <- dat$test_combo[which((dat$term == "knockout") )]
    } else {
        stop(paste(export, "is not a valid parameter"))
    }

    # order test-combos by estimate size for knockouts
    for_ordering <- dat %>% 
      filter(test_combo %in% combos_plot & term == "knockout") 

    # if one does not exists use the other
    bool_cond <- is.na(for_ordering$estimate_conditioned)
    for_ordering$estimate_conditioned[bool_cond] <-  for_ordering$estimate_unconditioned[bool_cond]


    for_ordering$p_to_use <- unlist(ifelse(!is.na(for_ordering$std.error_conditioned), 
                                           (for_ordering$p.value_conditioned), 
                                           (for_ordering$p.value_unconditioned)))
    for_ordering$est_plot <- for_ordering$p_to_use
    #for_ordering$est_plot <- pmax(for_ordering$p.value_unconditioned, 
    #                              for_ordering$p.value_conditioned)

    for_ordering <- for_ordering %>%
      arrange(est_plot, desc = F) %>%
      mutate(print_combo = paste0(hgnc_symbol, ":", diagnosis))


    # Only plot knockout HRs
    forest_dat <- bind_rows(uncond, cond) 
    forest_dat$hgnc_symbol <- mapping[forest_dat$gene]
    forest_dat <- forest_dat %>%
      filter(test_combo %in% combos_plot & term %in% c("wt","knockout","chet_cis")) %>%
      mutate(lci = estimate - 1.96*std.error,
             uci = estimate + 1.96*std.error,
             analysis_type = factor(analysis_type, levels = c("conditioned", "unconditioned")),
             plabel = paste0("P = ", gsub("e", "E", signif(p.value, digits = 3))),
             print_combo = factor(paste0(hgnc_symbol, ":", diagnosis), 
                                  levels = rev(unique(for_ordering$print_combo))))

    # for writing extended data tables
    outtable <- forest_dat
    outtable$vs <- "het"
    outtable <- outtable[,c("gene","hgnc_symbol", "term","vs", "diagnosis", "estimate","std.error", "p.value", "lci", "uci", "FDR", "analysis_type")]
    
    outfile1 <- paste0(out_prefix, ".txt.gz")
    write(paste("writing", outfile1), stdout())
    fwrite(outtable, outfile1, sep="\t")
    outfile2 <- paste0(out_prefix, ".ko.txt.gz")
    outtable <- outtable[outtable$term == "knockout",]
    write(paste("writing", outfile2), stdout())
    fwrite(outtable, outfile2, sep="\t")


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir_cond", default=NULL, required = TRUE, help = "")
parser$add_argument("--in_dir_uncond", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
parser$add_argument("--n_tests", default=NULL, required = TRUE, help = "")
parser$add_argument("--export", default=NULL, required = FALSE, help = "either NULL or 'everything'")
args <- parser$parse_args()

main(args)











