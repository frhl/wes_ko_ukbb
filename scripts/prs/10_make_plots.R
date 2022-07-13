library(ggplot2)
library(dplyr)
library(ggsci)
library(ggrastr)
library(ggrepel)
library(data.table)
library(latex2exp)
library(argparse)

main <- function(args){

    print(args)
    stopifnot(file.exists(args$summary_ldsc)) 
    stopifnot(file.exists(args$summary_prs_cts))
    stopifnot(file.exists(args$summary_prs_bin))
    stopifnot(dir.exists(dirname(args$out_prefix)))
 
    ldsc <- fread(args$summary_ldsc)
    cts <- fread(args$summary_prs_cts)
    bin <- fread(args$summary_prs_bin) 

    # read in relevant headers of phenotypes
    bin_header <- unlist(strsplit(args$bin_header, split = ","))
    cts_header <- unlist(strsplit(args$cts_header, split = ","))
    
    binary_header <- readLines(args$bin_header)
    cts_header <- paste0(readLines(args$cts_header),"_int")

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

    # ** cts traits **
    p1 <- ggplot(res_cts[res_cts$ldsc_coef == "h2",], 
       aes(
           y=phenotype, #reorder(phenotype, ldsc_estimate),
           x=ldsc_estimate, 
           xmax=ldsc_estimate+ldsc_std_error,
           xmin=ldsc_estimate-ldsc_std_error,
           color=included
           )
      ) + 
    geom_pointrange() +
    xlab(expression(~h^2 ~ "(SE)")) + 
    ylab("Phenotypes") +
    theme_bw() +
    scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
    scale_color_d3('category20c') +
    coord_cartesian(xlim=c(0, 0.6)) +
    guides(color="none")

    p2 <- ggplot(res_cts[res_cts$ldsc_coef == "intercept",], 
         aes(
             y=phenotype,
             x=ldsc_estimate, 
             xmax=ldsc_estimate+ldsc_std_error,
             xmin=ldsc_estimate-ldsc_std_error,
             color=included
             )
        ) + 
      theme_bw() +
      geom_pointrange() +
      xlab(expression("intercept (SE)")) + 
      geom_vline(xintercept = 1, linetype = 'dashed') +
      scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
      scale_color_d3('category20c') + 
      theme(axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) +
      guides(color="none")

    p3 <- ggplot(res_cts[res_cts$ldsc_coef == "h2",], 
         aes(
             y=phenotype, #y=reorder(phenotype, ldsc_estimate),
             x=-log10(ldsc_pvalue),
             xmax=-log10(ldsc_pvalue),
             xmin=-log10(ldsc_pvalue),
             color=included
             )
        ) + 
      geom_pointrange() +
      xlab(expression(log[10]*"(P-value)")) +
      geom_vline(xintercept = -log10(1e-5), linetype = 'dashed') +
      theme_bw() +
      scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
      scale_color_d3('category20c') +
      theme(axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) +
      guides(color="none")

    p4 <- ggplot(res_cts, 
         aes(
             y=phenotype, #y=reorder(phenotype, ldsc_estimate),
             x=corr,
             xmax=corr_ci_upper,
             xmin=corr_ci_lower,
             color=included
             )
        ) + 
      geom_point() +
      geom_pointrange() +
      theme_bw() +
      scale_color_d3('category20c') +
      geom_vline(xintercept = 0, linetype = 'dashed') +
      xlab(expression("Pearson correlation (95% CI)")) + 
      ylab("Phenotypes") +
      scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
      theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
      guides(color="none")
 
    # save resulting as a cobmined plot     
    out_p1 <- paste0(args$out_prefix, "_cts.png")
    out_d1 <- paste0(args$out_prefix, "_cts.txt.gz")
    fwrite(res_cts, out_d1, sep = "\t")
    png(out_p1)
    plt <- cowplot::plot_grid(p1, p2, p3, p4, rel_widths = (c(0.45, 0.2, 0.2, 0.2)), ncol = 4, nrow = 1)
    print(plt)
    graphics.off()

    p1 <- ggplot(res_cts[res_cts$ldsc_coef == "h2",], 
       aes(
           y=phenotype, #reorder(phenotype, ldsc_estimate),
           x=ldsc_estimate, 
           xmax=ldsc_estimate+ldsc_std_error,
           xmin=ldsc_estimate-ldsc_std_error,
           color=included
           )
        ) + 
        geom_pointrange() +
        xlab(expression(~h^2 ~ "(SE)")) + 
        ylab("Phenotypes") +
        theme_bw() +
        scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
        scale_color_d3('category20c') +
        coord_cartesian(xlim=c(0, 0.6)) +
        guides(color="none")

    p2 <- ggplot(res_cts[res_cts$ldsc_coef == "intercept",], 
           aes(
               y=phenotype,
               x=ldsc_estimate, 
               xmax=ldsc_estimate+ldsc_std_error,
               xmin=ldsc_estimate-ldsc_std_error,
               color=included
               )
          ) + 
        theme_bw() +
        geom_pointrange() +
        xlab(expression("intercept (SE)")) + 
        geom_vline(xintercept = 1, linetype = 'dashed') +
        scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
        scale_color_d3('category20c') + 
        theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
        guides(color="none")

    p3 <- ggplot(res_cts[res_cts$ldsc_coef == "h2",], 
           aes(
               y=phenotype, #y=reorder(phenotype, ldsc_estimate),
               x=-log10(ldsc_pvalue),
               xmax=-log10(ldsc_pvalue),
               xmin=-log10(ldsc_pvalue),
               color=included
               )
          ) + 
        geom_pointrange() +
        xlab(expression(log[10]*"(P-value)")) +
        geom_vline(xintercept = -log10(1e-5), linetype = 'dashed') +
        theme_bw() +
        scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
        scale_color_d3('category20c') +
        theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
        guides(color="none")

    p4 <- ggplot(res_cts, 
           aes(
               y=phenotype, #y=reorder(phenotype, ldsc_estimate),
               x=corr,
               xmax=corr_ci_upper,
               xmin=corr_ci_lower,
               color=included
               )
          ) + 
        geom_point() +
        geom_pointrange() +
        theme_bw() +
        scale_color_d3('category20c') +
        geom_vline(xintercept = 0, linetype = 'dashed') +
        xlab(expression("Pearson correlation (95% CI)")) + 
        ylab("Phenotypes") +
        scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
        theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
        guides(color="none")

    # save resulting as a cobmined plot     
    out_p2 <- paste0(args$out_prefix, "_bin.png")
    out_d2 <- paste0(args$out_prefix, "_bin.txt.gz")
    fwrite(res_bin, out_d1, sep = "\t")
    png(out_p2, width = as.numeric(args$out_width), height = as.numeric(args$out_height))
    plt <- cowplot::plot_grid(p1, p2, p3, p4, rel_widths = (c(0.45, 0.2, 0.2, 0.2)), ncol = 4, nrow = 1)
    print(plt)
    graphics.off()

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--summary_ldsc", default=NULL, required = TRUE, help = "LDSC summary file")
parser$add_argument("--summary_prs_cts", default=NULL, required = TRUE, help = "PRS cts validation file")
parser$add_argument("--summary_prs_bin", default=NULL, required = TRUE, help = "PRS binary validation file")
parser$add_argument("--bin_header", default=NULL, required = TRUE, help = "comma-seperated list of phenotypes")
parser$add_argument("--cts_header", default=NULL, required = TRUE, help = "comma-seperated list of phenotypes")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--out_width", default=12, required = FALSE, help = "plot size")
parser$add_argument("--out_height", default=14, required = FALSE, help = "plot size")
args <- parser$parse_args()

main(args)









