# define consequence
PLOF_CSQS = c("transcript_ablation", "splice_acceptor_variant","splice_donor_variant", "stop_gained", "frameshift_variant")
MISSENSE_CSQS = c("stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant")
SYNONYMOUS_CSQS = c("stop_retained_variant", "synonymous_variant")
OTHER_CSQS = c("mature_miRNA_variant", "5_prime_UTR_variant","3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant","NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant","downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant","regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation","regulatory_region_variant", "feature_truncation", "intergenic_variant")
                                
# add consequence category
annotate_csqs_category <- function(csqs){
    csqs[csqs %in% PLOF_CSQS] <- "ptv"
    csqs[csqs %in% MISSENSE_CSQS] <- "missense"
    csqs[csqs %in% SYNONYMOUS_CSQS] <- "synonymous"
    csqs[csqs %in% OTHER_CSQS] <- "non_coding"
    return(csqs)
}     

