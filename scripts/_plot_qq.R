#!/usr/bin/env Rscript

library(argparse)
library(ggsci)
library(ggplot2)
library(reshape2)
library(data.table)

# remove NULLs from list
null_omit <- function(lst) lst[-which(sapply(lst, is.null))]

main <- function(args){

  stopifnot(file.exists(args$bridge))
  stopifnot(file.exists(args$phenotypes))

  # plotting parameters
  n_x_ticks=10 
  n_y_ticks=10
  alpha = 0.9 

  # Load bridge between hgnc_symbol and ensembl
  bridge <- fread(args$bridge)
  ensembl_to_hgnc <- bridge$hgnc_symbol
  ensembl_to_chr <- bridge$chromosome_name
  ensembl_to_pos <- round((bridge$start_position + bridge$end_position) / 2)
  names(ensembl_to_hgnc) <- bridge$ensembl_gene_id
  names(ensembl_to_pos) <- bridge$ensembl_gene_id
  names(ensembl_to_chr) <- bridge$ensembl_gene_id

  # load phenotypes
  phenotypes <- fread(args$phenotypes)
  phenos <- colnames(phenotypes)

  lst <- lapply(phenos, function(ph) {
    
    # retrieve path for SAIGE analysis
    #path_normal <- paste0("data/saige/output/",trait,"/step2/min_mac4/ukb_eur_wes_200k_maf0to5e-2_",ph,"_pLoF_damaging_missense.txt.gz")
    #path_prs <- paste0("data/saige/output/",trait,"/step2/min_mac4/ukb_eur_wes_200k_maf0to5e-2_",ph,"_pLoF_damaging_missense_locoprs.txt.gz" ) 
    
    # Substiute string to get the right phenotype
    path_normal <- gsub("PHENOTYPE", ph, args$path_marker_normal)
    path_prs <- gsub("PHENOTYPE", ph, args$path_marker_prs)    

    # check paths
    path_normal_exists <- file.exists(path_normal)
    path_prs_exists <- file.exists(path_prs)
    
    # read single marker file with/without PRS
    if (path_normal_exists){
        cols <- c("MarkerID","p.value")
        d_normal <- fread(path_normal)
        d_normal <- d_normal[,colnames(d_normal) %in% cols, with = FALSE]
        # plot PRS as well if present
        if (path_prs_exists){
            d_prs <- fread(path_prs)
            d_prs <- d_prs[,colnames(d_prs) %in% cols, with = FALSE]
            colnames(d_prs)[2] <- "p.value_prs"
            d <- merge(d_normal, d_prs, by = "MarkerID")
            colnames(d)[2:3] <- c("without PRS","with PRS")
            dt <- melt(d, 
                       id.vars = c("MarkerID"), 
                       variable.name = "test", 
                       value.name = "Pvalue")
        } else {
            dt <- d_normal
            colnames(dt) <- c("MarkerID", "Pvalue")
            dt$test <- "without PRS"   
        }
        # obtain 95% confidence intervals for qq-plots
        dt$Pvalue <- -log10(dt$Pvalue)
        dt <- data.table(
            dt %>% arrange(test, desc(Pvalue)) %>% select(-test),
            dt %>% group_by(test) %>% arrange(desc(Pvalue)) %>% 
                summarize(
                Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
                clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
                cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
                )
            )
        # how many markers labels should be drawn
        markers <- head(dt[rev(order(dt$Pvalue))]$MarkerID, n = 3)
        # append some info to data.table
        dt$phenotype <- ph
        dt$hgnc_symbol <- ensembl_to_hgnc[dt$MarkerID]
        dt$chrom <- ensembl_to_chr[dt$MarkerID]
        dt$position <- ensembl_to_pos[dt$MarkerID]
        dt$label <- NA
        dt$label[dt$MarkerID %in% markers] <- dt$hgnc_symbol[dt$MarkerID %in% markers]
        dt$label_phenotype <- paste0(dt$hgnc_symbol,'-',dt$phenotype)
        dt$label_phenotype[is.na(dt$label)] <- NA
        return(dt)
    }
  })

  # combine all data.points
  lst <- null_omit(lst)
  n_plots <- length(lst)
  dts <- do.call(rbind, lst)

  p <- ggplot(dts, aes(x=Pvalue_expected, y = Pvalue, ymin=clower, ymax=cupper, color = test, label = label)) +
    geom_ribbon(fill="grey80", color="grey80") +
    geom_hline(yintercept=-log10(p_threshold), linetype = 'dashed', color = 'red') +
    geom_point_rast(alpha=alpha, raster.dpi=500) +
    geom_abline(linetype = 'dashed') +
    geom_label_repel(
        box.padding = 0.4, label.padding=0.25, point.padding = 0.2,
        color = 'grey30', segment.color = 'grey50', max.overlaps=Inf,
        size=4, segment.size=0.1, show.legend = FALSE
    ) + 
    theme_bw() +
    xlab(expression(-log[10]*"(P-value expected)")) +
    ylab(expression(-log[10]*"(P-value observed)")) +
    scale_color_d3('category20c', limits=NULL) +
    labs(color=paste0("Analysis")) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=n_x_ticks)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=n_y_ticks)) +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
          axis.title.y = element_text(margin=ggplot2::margin(r=10)),
          plot.title = element_text(hjust=0.5),
          legend.position="bottom") 

   # setup output plotting
   facet_rows <- args$facet_rows
   facet_cols <- args$facet_cols 
   plots_per_page <- facet_rows * facet_col
   pages <- ceiling(n_plots / plots_per_pags)

   # plot results
   outpdf <- paste0(args$prefix,'.pdf')
   pdf(outpdf, width = 12, height = 16)
   for (p in 1:pages){
    plt <- p + facet_wrap_paginate(~phenotype, scales = "free", nrow = facet_rows, ncol = facet_cols, page = p)
    out_page <- paste0(args$prefix, "_",p,"of",pages,".png")
    ggsave(out_page, plt, width = 12, height = 16)
    print(plt)
   } 
   graphics.off()

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--p_threshold", default=5e-6, help = "?")
parser$add_argument("--bridge", default="/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz", help = "?")
parser$add_argument("--phenotypes", default=NULL, help = "?")
parser$add_argument("--path_marker_normal", default=NULL, help = "?")
parser$add_argument("--path_marker_prs", default=NULL, help = "?")
parser$add_argument("--facet_rows", default=5, help = "?")
parser$add_argument("--facet_cols", default=4, help = "?")
args <- parser$parse_args()

main(args)

