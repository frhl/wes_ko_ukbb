#!/usr/bin/env Rscript

library(argparse)
library(data.table)


plof_csqs = c("transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant")

missense_csqs = c("stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant")

synonymous_csqs = c("stop_retained_variant", "synonymous_variant")

other_csqs = c("mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant")


# get counts of variants by category
tabulate_variants <- function(chr, vep_path, vep_min_mac_path = NULL){
    # read in data
    path_chr <- gsub("CHR", chr, vep_path)
    d <- fread(path_chr)
    if (!is.null(vep_min_mac_path)){
        path_vep_min_mac <- gsub("CHR", chr, vep_min_mac_path)
        d_vep_min_mac <- fread(path_vep_min_mac)
        # create dict for min mac mapping
        stopifnot(d_vep_min_mac$rsid %in% d$rsid)
        min_mac_dict <- d_vep_min_mac$MAC
        names(min_mac_dict) <- d_vep_min_mac$rsid
        # map min mac
        d$MAC_old <- d$MAC
        d$MAC <- min_mac_dict[d$rsid]
        # subset so that no rows with zero variants are kept
        
    }
    d <- d[d$MAC>0,]
    d$is_singleton <- d$MAC == 1
    df_all <- as.data.frame(
        table(
            d$csqs.most_severe_consequence,
            d$consequence_category
        )
    )
    df_singletons <- as.data.frame(
        table(
            d$csqs.most_severe_consequence[d$is_singleton],
            d$consequence_category[d$is_singleton]
        )
    )
    # setup column names
    colnames(df_all)[3] <- paste0("n_total_chr", chr)
    colnames(df_singletons)[3] <- paste0("n_singletons_chr", chr)
    # combine singletons and total
    df <- setDT(merge(df_all, df_singletons, by = c("Var1", "Var2"), all = TRUE))
    colnames(df)[1] <- "variant_category"
    colnames(df)[2] <- "consequence_category"
    return(df)
    
}

# read in all chromosomes and combine counts of variants by category
get_vep_table <- function(vep_path, vep_min_mac_path = NULL){
    if (!is.null(vep_min_mac_path)) write(paste("Using min_mac from ", basename(vep_min_mac_path)), stdout())
    autosomes <- 1:22
    lst_variant_category <- lapply(autosomes, function(chr){
      return(tabulate_variants(chr, vep_path, vep_min_mac_path))
    })
    # write table of variant consequences
    dt_variant_category <- Reduce(merge, lst_variant_category)
    dt_variant_category[is.na(dt_variant_category)] <- 0
    return(dt_variant_category)
}

# create super categories
create_variant_table <- function(dt_variant_category){
    # combine into single table
    dt_aggr_variant_category <- data.table(
        variant_category = dt_variant_category$variant_category,
        consequence_category = dt_variant_category$consequence_category,
        n_total = rowSums(dt_variant_category[,grepl("n_total",colnames(dt_variant_category)),with=FALSE]),
        n_singletons = rowSums(dt_variant_category[,grepl("n_singletons",colnames(dt_variant_category)),with=FALSE])
    )
    # count fraction of singletons
    dt_aggr_variant_category <- dt_aggr_variant_category[dt_aggr_variant_category$n_total > 0,]
    dt_aggr_variant_category$pct_singletons <- round((dt_aggr_variant_category$n_singletons / dt_aggr_variant_category$n_total)*100, 1)
    # create super category
    dt_aggr_variant_category$super_category <- NA
    dt_aggr_variant_category$super_category[dt_aggr_variant_category$variant_category %in% plof_csqs] <- "PTV"
    dt_aggr_variant_category$super_category[dt_aggr_variant_category$variant_category %in% missense_csqs] <- "Missense"
    dt_aggr_variant_category$super_category[dt_aggr_variant_category$variant_category %in% synonymous_csqs] <- "Synonymous"
    dt_aggr_variant_category$super_category[dt_aggr_variant_category$variant_category %in% other_csqs] <- "Non-coding"
    dt_aggr_variant_category$super_category[is.na(dt_aggr_variant_category$super_category)] <- "Non-coding"
    dt_aggr_variant_category = dt_aggr_variant_category[order(dt_aggr_variant_category$super_category)]
    dt_aggr_variant_category <- dt_aggr_variant_category[,c(6,2,1,3,4,5)]
    dt <- dt_aggr_variant_category
    return(dt)
}

# clean up and make pretty
make_variant_table_pretty <- function(dt){
    categories <- c("PTV", "Missense","Synonymous", "Non-coding")
    lst <- lapply(categories, function(category){
        dt_subset <- dt[dt$super_category %in% category,]
        # get combined data
        dt_total_cat <- aggregate(n_total ~ consequence_category, data = dt_subset, FUN = sum)
        dt_singleton_cat <- aggregate(n_singletons ~ consequence_category, data = dt_subset, FUN = sum)
        # aggregate and merge
        dt_cat_mrg <- merge(dt_total_cat, dt_singleton_cat)
        dt_cat_mrg$variant_category <- " "
        dt_cat_mrg$super_category <- category
        dt_cat_mrg$pct_singletons <- round((dt_cat_mrg$n_singletons / dt_cat_mrg$n_total)*100, 1)
        dt_cat_mrg <- dt_cat_mrg[,c(5,1,4,2,3,6)]
        # create total counts
        dt_total <- data.table(
            super_category = category,
            consequence_category = " ",
            variant_category = " ",
            n_total = sum(dt_cat_mrg$n_total),
            n_singletons = sum(dt_cat_mrg$n_singletons),
            pct_singletons = round((sum(dt_cat_mrg$n_singletons) / sum(dt_cat_mrg$n_total))*100, 1)
        )
        # combine with subset
        dt_out <- rbind(dt_subset, dt_cat_mrg, dt_total)
        return(dt_out)
    })
    outdt <- do.call(rbind,lst)
    return(outdt)
}




main <- function(args){

    in_file_prefilter <- args$in_file_prefilter
    in_file_postfilter <- args$in_file_postfilter
    out_prefix <- args$out_prefix

    # pre-filtered VEP variants
    dt1 <- get_vep_table(in_file_prefilter)
    dt1 <- create_variant_table(dt1)
    dt1 <- make_variant_table_pretty(dt1)
    outfile1 <- paste0(out_prefix, "_prefilter.txt.gz")
    fwrite(dt1, outfile1, sep = "\t")
     
    # pre and post-filtered VEP variants (note that we use 
    # the MAC from the old file to figure out exactly how 
    # many of previous singletons are retained)
    dt2 <- get_vep_table(in_file_postfilter, in_file_prefilter)
    dt2 <- create_variant_table(dt2)
    dt2 <- make_variant_table_pretty(dt2)
    outfile2 <- paste0(out_prefix, "_postfilter.txt.gz")
    fwrite(dt2, outfile2, sep = "\t")
     
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_file_prefilter", default=NULL, help = "infile")
parser$add_argument("--in_file_postfilter", default=NULL, help = "infile")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

