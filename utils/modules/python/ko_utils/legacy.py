#!/usr/bin/env python3

def annotate_vep(mt, vep_path):
    r'''Annotate matrix table with VEP consequence from external file.'''
    print(f'Annotating with VEP file: {vep_path}')
    
    # Open file containing VEP fields
    with open('data/vep/vep_fields.txt', 'r') as file:
        fields = file.read().strip().split(',')
    ht = hl.import_vcf(vep_path).rename({'info':'vep'}) 
    
    # Add VEP fields by iteration
    for i in range(len(fields)):
        ht = ht.annotate_rows(
            vep=ht.vep.annotate(
                col=ht.vep.CSQ.map(lambda x: (x.split('\\|')[i]))[0]
                ).rename({'col':f'{fields[i]}'})
        )
    
    # Most severe variant consequence
    ht = ht.annotate_rows(vep = ht.vep.annotate(most_severe_consequence = ht.vep.Consequence.split('&')[0]))
    
    # Extract various categories annotations and change type
    ht = ht.annotate_rows(vep = ht.vep.annotate(sift_pred = ht.vep.SIFT_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hdiv_pred = ht.vep.Polyphen2_HDIV_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(polyphen2_hvar_pred = ht.vep.Polyphen2_HVAR_pred.split('&')[0]))
    ht = ht.annotate_rows(vep = ht.vep.annotate(cadd_phred_score = hl.parse_float(ht.vep.CADD_phred)))
    ht = ht.annotate_rows(vep = ht.vep.annotate(revel_score = hl.parse_float(ht.vep.REVEL_score)))
    
    # Define protein truncating variants
    ptv = hl.set(["transcript_ablation", "splice_acceptor_variant",
              "splice_donor_variant", "stop_gained", "frameshift_variant"])
    
    # Define missense variation
    missense = hl.set(["stop_lost", "start_lost", "transcript_amplification",
                   "inframe_insertion", "inframe_deletion", "missense_variant",
                   "protein_altering_variant", "splice_region_variant"])
    
    # Define synonymous
    synonymous = hl.set(["incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant"])
    
    # Define non coding variation
    non_coding = hl.set(["coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant"])
    
    # Create categories for downstream analysis
    ht = ht.annotate_rows(vep = ht.vep.annotate(consequence_category = 
        hl.case().when(ptv.contains(ht.vep.most_severe_consequence), "ptv")
             .when(missense.contains(ht.vep.most_severe_consequence) & 
                   (~hl.is_defined(ht.vep.cadd_phred_score) | 
                    ~hl.is_defined(ht.vep.revel_score)), "other_missense")                                   
             .when(missense.contains(ht.vep.most_severe_consequence) & 
                   (ht.vep.cadd_phred_score >= 20) & 
                   (ht.vep.revel_score >= 0.6), "damaging_missense") 
             .when(missense.contains(ht.vep.most_severe_consequence), "other_missense")
             .when(synonymous.contains(ht.vep.most_severe_consequence), "synonymous")
             .when(non_coding.contains(ht.vep.most_severe_consequence), "non_coding")
             .default("NA")
    ))
                                                
    # combine with matrix table
    mt = mt.annotate_rows(vep = ht.index_rows(mt.locus, mt.alleles).vep.drop('CSQ'))
    return(mt)

