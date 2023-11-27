#!/usr/bin/env python3

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import io
from ko_utils import ko

def main(args):

    # parser
    vep_path = args.vep_path
    out_prefix = args.out_prefix
    spliceai_path = args.spliceai_path
    spliceai_score = float(args.spliceai_score)
    revel_score = float(args.revel_score)
    cadd_score = float(args.cadd_score)  # Fixed: use cadd_score instead of revel_score
    case_builder = args.case_builder 
    
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    ht = hl.read_table(vep_path)

    # key variant (snp, coding sequence id)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(transcript_id=ht.vep.worst_csq_by_gene_canonical.transcript_id)
    #ht = ht.annotate(ccds=ht.vep.worst_csq_by_gene_canonical.ccds.split("\\.")[0])
    #ht = ht.key_by(ht.locus, ht.alleles, ht.ccds) 
    ht = ht.key_by(ht.locus, ht.alleles, ht.transcript_id) 

    # get alpha missense and key variant by snp and ccds
    #alpha = hl.read_table(alpha_missense_path)
    #alpha = alpha.key_by(alpha.locus, alpha.alleles, alpha.ccds)
    #alpha = alpha.semi_join(ht)

    #ht = ht.annotate(alpha_info=alpha[ht.key])
    # annotate consequence
    #ht = ht.transmute(
    #    vep = ht.vep.annotate(
    #        worst_csq_by_gene_canonical = ht.vep.worst_csq_by_gene_canonical.annotate(
    #            am_pathogenicity = hl.float(ht.alpha_info.am_pathogenicity),
    #            am_class = ht.alpha_info.am_class
    #        )
    #    )
    #)

    # get spliceAI annotation, note that we if we don't key by
    # mane select, then we get downstream gene variants annotated as splice variantss
    spliceai = hl.read_table(spliceai_path)
    spliceai = spliceai.annotate(
            transcript_id=spliceai.SpliceAI.SYMBOL.split("---")[2].split("\\.")[0]
    )
   
    # only keep overlapping entries
    spliceai = spliceai.key_by(spliceai.locus, spliceai.alleles, spliceai.transcript_id)
    spliceai = spliceai.semi_join(ht)

    # index into spliceAI and extract annotaions     
    ht = ht.annotate(spliceai_info=spliceai[ht.key])
    ht = ht.transmute(
        vep = ht.vep.annotate(
            worst_csq_by_gene_canonical = ht.vep.worst_csq_by_gene_canonical.annotate(
                SpliceAI_DS_max = ht.spliceai_info.SpliceAI.DS_max
            )
        )
    )

    print(ht.count())
    # drop the MANE columns
    ht = ht.key_by(ht.locus, ht.alleles)
    ht = ht.drop(ht.transcript_id)

    #if case_builder in "alpha_missense":
    #    ht = ht.annotate(
    #        brava_csqs=ko.csqs_case_builder_alpha_missense_v2(
    #                worst_csq_expr=ht.vep.worst_csq_by_gene_canonical,
    #                spliceai_cutoff = spliceai_score,
    #                revel_cutoff = revel_score,
    #                am_cutoff = am_score
    #        )
    #    )
    if case_builder in "brava":
        print(f"Using {case_builder} builder with cadd>{cadd_score}\trevel>{revel_score}\tspliceAI>{spliceai_score}")
        ht = ht.annotate(
                    brava_csqs=ko.csqs_case_builder_brava(
                            worst_csq_expr=ht.vep.worst_csq_by_gene_canonical,
                            spliceai_cutoff = spliceai_score,
                            revel_cutoff = revel_score,
                            cadd_cutoff = cadd_score
             )
         )
    elif case_builder in "original":
         ht = ht.annotate(
                    brava_csqs=ko.csqs_case_builder(
                            worst_csq_expr=ht.vep.worst_csq_by_gene_canonical
             )
         )
    else:
        raise ValueError("Not a valid case_builder argument!")


    # quick annotated with some useful info
    ht = ht.annotate(
            gene_symbol=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
            gene_id=ht.vep.worst_csq_by_gene_canonical.gene_id,
            transcript=ht.vep.worst_csq_by_gene_canonical.transcript_id,
            biotype=ht.vep.worst_csq_by_gene_canonical.biotype,
            mane_select=ht.vep.worst_csq_by_gene_canonical.mane_select,
            canonical=ht.vep.worst_csq_by_gene_canonical.canonical,
            ccds=ht.vep.worst_csq_by_gene_canonical.ccds,
            csqs=ht.vep.worst_csq_by_gene_canonical.most_severe_consequence,
            revel_score=ht.vep.worst_csq_by_gene_canonical.revel_score,
            cadd_phred=ht.vep.worst_csq_by_gene_canonical.cadd_phred,
            spliceai_max_ds=ht.vep.worst_csq_by_gene_canonical.SpliceAI_DS_max,
            loftee_lof=ht.vep.worst_csq_by_gene_canonical.lof
            #am_class=ht.vep.worst_csq_by_gene_canonical.am_class,
            #am_pathogenicity=ht.vep.worst_csq_by_gene_canonical.am_pathogenicity
    )

    # filter to protein coding
    ht = ht.filter(
        ht.biotype == "protein_coding"
    )

    # annotate with actual variant ID
    ht = ht.annotate(
        varid=hl.delimit([
            hl.str(ht.locus.contig),
            hl.str(ht.locus.position),
            ht.alleles[0],
            ht.alleles[1]],':')
    )

    ht = ht.select(*[ht.varid, ht.gene_symbol, ht.gene_id, ht.transcript, ht.biotype, ht.mane_select, ht.canonical, ht.csqs, ht.brava_csqs, ht.revel_score,ht.cadd_phred,ht.loftee_lof,ht.spliceai_max_ds])
    ht.write(out_prefix + ".ht", overwrite=True)
    ht.export(out_prefix + ".txt.gz")


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--case_builder', default=None, help='What case builder should be used')
    parser.add_argument('--vep_path', default=None, help='Path to input')
    parser.add_argument('--alpha_missense_path', default=None, help='Path to spliceai VCF')
    parser.add_argument('--spliceai_path', default=None, help='Path to spliceai VCF')
    parser.add_argument('--spliceai_score', type=float, default=0.50, help='SpliceAI delta score')
    parser.add_argument('--revel_score', type=float, default=0.773, help='Revel score cutoff')
    parser.add_argument('--cadd_score', type=float, default=28.1, help='CADD score cutoff')
    #parser.add_argument('--am_score', type=float, default=0.564, help='AlphaMissense score cutoff')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)

