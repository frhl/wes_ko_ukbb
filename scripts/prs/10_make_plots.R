
library(argparse)
library(data.table)
library(bigstatsr)

main <- function(args){

    print(args)
    stopifnot(file.exists(args$summary_ldsc)) 
    stopifnot(file.exists(args$summary_prs_cts))
    stopifnot(file.exists(args$summary_prs_bin))
    stopifnot(dir.exists(dirname(args$out_prefix)))
 
    ldsc <- fread(args$summary_ldsc)
    cts <- fread(args$summary_prs_cts)
    bin <- fread(args$summary_prs_bin) 

    bin_header <- unlist(strsplit(args$bin_header, split = ",")
    cts_header <- unlist(strsplit(args$cts_header, split = ",")
    
    # rename ldsc so that it can be used for merging
    colnames(ldsc) <- paste0("ldsc_",colnames(ldsc))
    colnames(ldsc)[colnames(ldsc) == "ldsc_phenotype"] <- "phenotype"

    # merge files
    res_bin <- merge(ldsc, bin, all.x = TRUE)
    res_bin <- res_bin[res_bin$phenotype %in% binary_header,]
    res_cts <- merge(ldsc, cts, all.x = TRUE)
    res_cts <- res_cts[res_cts$phenotype %in% cts_header,]

    # are they used for prs?
    res_bin$included <- ifelse(!is.na(res_bin$auc_mean), "Yes", "No")
    res_cts$included <- ifelse(!is.na(res_cts$correlation), "Yes", "No")

    stopifnot(nrow(res_cts) > 0)
    stopifnot(nrow(res_bin) > 0)

    # ** binary traits **
    p1 <- ggplot(res_bin[res_bin$ldsc_coef == "h2",], 
       aes(
           x=phenotype, #reorder(phenotype, ldsc_estimate),
           y=ldsc_estimate, 
           ymax=ldsc_estimate+ldsc_std_error,
           ymin=ldsc_estimate-ldsc_std_error,
           fill = included
           )
      ) + 
      geom_bar(stat="identity", position = 'dodge') +
      geom_point() +
      geom_errorbar() +
      ylab(bquote(~h^2)) + 
      xlab("") +
      labs(fill="PRS") +
      theme_bw() +
      coord_cartesian(ylim=c(0, 1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    p2 <- ggplot(res_bin[!is.na(res_bin$auc_mean),], 
       aes(
           x=phenotype, 
           y=auc_mean,
           ymax=auc_97_5_pct,
           ymin=auc_2_5_pct
           )
      ) + 
      geom_point() +
      geom_errorbar() +
      geom_hline(yintercept = 0.5, linetype = 'dashed') +
      xlab("") + 
      ylab("AUC") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    out_p1 <- paste0(args$out_prefix,"_bin_h2.png")
    out_p2 <- paste0(args$out_prefix,"_bin_auc.png")
    ggsave(p1, out_p1, width = 14, height = 7)
    ggsave(p2, out_p2, width = 7, height = 6)

    p3 <- ggplot(res_cts[res_cts$ldsc_coef == "h2",], 
           aes(
               x=phenotype, #reorder(phenotype, ldsc_estimate),
               y=ldsc_estimate, 
               ymax=ldsc_estimate+ldsc_std_error,
               ymin=ldsc_estimate-ldsc_std_error,
               fill = included
               )
          ) + 
        geom_bar(stat="identity", position = 'dodge') +
        geom_point() +
        geom_errorbar() +
        ylab(bquote(~h^2)) + 
        xlab("") +
        labs(fill="PRS") +
        theme_bw() +
        coord_cartesian(ylim=c(0, 1)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    p4 <- ggplot(res_cts[!is.na(res_cts$correlation),], 
           aes(
               x=phenotype, 
               y=correlation
               )
          ) + 
        geom_point() +
        geom_hline(yintercept = 0, linetype = 'dashed') +
        xlab("") + 
        ylab("Pearson R2") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
 
    out_p3 <- paste0(args$out_prefix,"_cts_h2.png")
    out_p4 <- paste0(args$out_prefix,"_cts_pearson.png")
    ggsave(p3, out_p3, width = 14, height = 7)
    ggsave(p4, out_p4, width = 7, height = 6)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--summary_ldsc", default=NULL, required = TRUE, help = "LDSC summary file")
parser$add_argument("--summary_prs_cts", default=NULL, required = TRUE, help = "PRS cts validation file")
parser$add_argument("--summary_prs_bin", default=NULL, required = TRUE, help = "PRS binary validation file")
parser$add_argument("--bin_header", default=NULL, required = TRUE, help = "comma-seperated list of phenotypes")
parser$add_argument("--cts_header", default=NULL, required = TRUE, help = "comma-seperated list of phenotypes")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









