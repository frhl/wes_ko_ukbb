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
    
    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
    ht = hl.read_table(vep_path)

    # key variant (snp, coding sequence id)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(transcript_id=ht.vep.worst_csq_by_gene_canonical.transcript_id)
    ht = ht.key_by(ht.locus, ht.alleles, ht.transcript_id) 

    # drop the MANE columns
    ht = ht.key_by(ht.locus, ht.alleles)
    ht = ht.drop(ht.transcript_id)

    ht = ht.annotate(
                brava_csqs=ko.csqs_case_builder(
                        worst_csq_expr=ht.vep.worst_csq_by_gene_canonical
         )
     )

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
            loftee_lof=ht.vep.worst_csq_by_gene_canonical.lof
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

    ht = ht.select(*[ht.varid, ht.gene_symbol, ht.gene_id, ht.transcript, ht.biotype, ht.mane_select, ht.canonical, ht.csqs, ht.brava_csqs, ht.revel_score,ht.cadd_phred,ht.loftee_lof])
    ht.write(out_prefix + ".ht", overwrite=True)
    ht.export(out_prefix + ".txt.gz")


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    # initial params
    parser.add_argument('--vep_path', default=None, help='Path to input')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')

    args = parser.parse_args()

    main(args)

