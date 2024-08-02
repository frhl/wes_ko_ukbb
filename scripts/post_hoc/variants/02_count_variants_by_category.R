#!/usr/bin/env Rscript

library(argparse)
library(ggsci)
library(ggplot2)
library(reshape2)
library(data.table)

super_categories <- list(
  PTV = c("transcript_ablation", "splice_acceptor_variant",
          "splice_donor_variant", "stop_gained", "frameshift_variant"),
  Missense = c("stop_lost", "start_lost", "transcript_amplification",
               "inframe_insertion", "inframe_deletion", "missense_variant"),
  Synonymous = c("stop_retained_variant", "synonymous_variant"),
  Non_coding = c("mature_miRNA_variant", "5_prime_UTR_variant",
                 "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
                 "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
                 "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
                 "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
                 "regulatory_region_variant", "feature_truncation", "intergenic_variant")
)

# Helper function to determine the super category of a variant
get_super_category <- function(csq) {
  for (category_name in names(super_categories)) {
    if (csq %in% super_categories[[category_name]]) {
      return(category_name)
    }
  }
  return('Non_coding')  # Default category if none match
}

get_total_row <- function(combined, cat_text, pred_text){
    
    combined <- combined[combined$summary_row==FALSE,]
    sum_n_variants_before_pp <- sum(combined$n_variants_before_pp)
    sum_n_singletons_before_pp <- sum(combined$n_singletons_before_pp)
    summed_pct_singletons_before_pp <- round((sum_n_singletons_before_pp / sum_n_variants_before_pp)*100, 1)

    sum_n_variants_after_pp <- sum(combined$n_variants_after_pp)
    sum_n_singletons_after_pp <- sum(combined$n_singletons_after_pp)
    summed_pct_singletons_after_pp <- round((sum_n_singletons_after_pp / sum_n_variants_after_pp)*100, 1)

    out <- data.table(
        category=cat_text,
        predicted_csqs=pred_text,
        vep_annotation="",
        n_variants_before_pp=sum_n_variants_before_pp,
        n_singletons_before_pp=sum_n_singletons_before_pp,
        pct_singletons_before_pp=summed_pct_singletons_before_pp,
        n_variants_after_pp=sum_n_variants_after_pp,
        n_singletons_after_pp=sum_n_singletons_after_pp,
        pct_singletons_after_pp=summed_pct_singletons_after_pp,
        summary_row=TRUE
    )  
    return(out)
}




main <- function(args){

    out_prefix <- args$out_prefix
    in_file <- args$in_file
    filter_biotype <- args$filter_biotype

    dt <- fread(in_file)
    if (!is.null(filter_biotype)) dt <- dt[dt$biotype == filter_biotype,]
    stopifnot(nrow(dt)>0)
    
    # add a few important covariates
    dt[, super_category := sapply(CSQ, get_super_category)]
    dt$is_singleton_before_pp <- dt$MAC.before_pp == 1
    dt$is_singleton_after_pp <- (dt$MAC.after_pp == 1) & (dt$is_singleton_before_pp)

    # get major categories and subset. This removes a few (non_coding CSQs)
    all_csqs <- as.character(unlist(super_categories))
    all_categories <- names(super_categories)
    dt <- dt[dt$super_category %in% all_categories,]
    dt <- dt[dt$CSQ %in% all_csqs,]  
    all_predictions <- unique(dt$annotation)

    out <- rbindlist(lapply(all_categories, function(cat){
        dt_by_cat <- dt[dt$super_category %in% cat, ]
        
        # middle loop for csqs prediction
        pred_loop <- rbindlist(lapply(all_predictions, function(pred){
            dt_by_cat_and_pred <- dt_by_cat[dt_by_cat$annotation %in% pred, ]
            
            # inner loop for VEP csqs
            csqs_loop <- rbindlist(lapply(all_csqs, function(csq){
                
                dt_subset_before_pp <- dt_by_cat_and_pred[dt_by_cat_and_pred$CSQ %in% csq,]
                n_variants_before_pp <- nrow(dt_subset_before_pp)
                n_singletons_before_pp <- sum(dt_subset_before_pp$is_singleton_before_pp)
                pct_singletons_before_pp <- round((n_singletons_before_pp / n_variants_before_pp)*100, 1)

                dt_subset_after_pp <- dt_subset_before_pp[dt_subset_before_pp$MAC.after_pp > 0,]
                n_variants_after_pp <- nrow(dt_subset_after_pp)
                n_singletons_after_pp <- sum(dt_subset_after_pp$is_singleton_before_pp)
                pct_singletons_after_pp <- round((n_singletons_after_pp / n_variants_after_pp)*100, 1)
                
                out <- data.table(
                    category=cat,
                    predicted_csqs=pred,
                    vep_annotation=csq,
                    n_variants_before_pp,
                    n_singletons_before_pp,
                    pct_singletons_before_pp,
                    n_variants_after_pp,
                    n_singletons_after_pp,
                    pct_singletons_after_pp,
                    summary_row=FALSE
                )
                return(out)
            }))
            
            # get total for the prediction 
            return(rbind(csqs_loop, get_total_row(csqs_loop, cat, paste("Total",pred))))
                    
        }))

        # get total for the full category
        return(rbind(pred_loop, get_total_row(pred_loop, paste("Total",cat), "")))
        
    }))

    out <- out[out$n_variants_before_pp>0,]
    outfile <- paste0(out_prefix, ".txt")
    write(paste("writing to", outfile), stdout())
    fwrite(out, outfile, sep="\t", quote=FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--filter_biotype", default=NULL, help = "?")
parser$add_argument("--in_file", default=NULL, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
args <- parser$parse_args()

main(args)


