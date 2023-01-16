#!/usr/bin/env python3

import hail as hl
import argparse

from ko_utils import io
from ukb_utils import hail_init
from ukb_utils import genotypes
from ukb_utils import variants
from ukb_utils import samples


def main(args):

    # parser
    hapmap = args.hapmap
    out_prefix = args.out_prefix

    hail_init.hail_bmrc_init_local('logs/hail/hail_format.log', 'GRCh38')
    hl._set_flags(no_whole_stage_codegen='1') # from zulip
   
    # load data and annotate
    ht = hl.import_table(hapmap, force = True, impute = True)
    ht = ht.annotate(grch37_varid = hl.delimit([hl.str(ht.chr), hl.str(ht.pos), hl.str(ht.a0), hl.str(ht.a1)], ':'))
    ht = ht.key_by(**hl.parse_variant(ht.grch37_varid, reference_genome = 'GRCh37'))

    # setup liftover references
    from_build = "GRCh37"
    to_build = "GRCh38"
    liftover_path = variants.get_liftover_chain_path(from_build, to_build)
    rg_from = hl.get_reference(from_build)
    rg_to = hl.get_reference(to_build)
    if not rg_from.has_liftover(rg_to):
        rg_from.add_liftover(liftover_path, rg_to)

    if not rg_to.has_sequence():
        reference_path = variants.get_reference_path("GRCh38")
        rg_to.add_sequence(*reference_path)

    # do liftover and flip allele for negative strand
    ht = ht.annotate(new_locus=hl.liftover(ht.locus, to_build, include_strand=True),old_locus=ht.locus)
    ht = ht.filter(hl.is_defined(ht.new_locus))
    ht = ht.annotate(new_alleles=hl.if_else(ht.new_locus.is_negative_strand, [hl.reverse_complement(ht.alleles[0]), hl.reverse_complement(ht.alleles[1])],ht.alleles))
    ht = ht.key_by(locus=ht.new_locus.result, alleles=ht.new_alleles)

    # check if the liftover was succesfull
    ht = ht.annotate(
        ref_allele_mismatch = ht.new_locus.result.sequence_context() != ht.new_alleles[0],
        locus_fail_liftover=hl.is_missing(ht.new_locus.result)
    )

    # annotate new locus
    ht = ht.annotate(grch38_varid = hl.delimit([hl.str(ht.locus.contig), hl.str(ht.locus.position), hl.str(ht.alleles[0]), hl.str(ht.alleles[1])], ':'))

    # export file
    ht = ht.select(*[ht.grch37_varid, ht.grch38_varid, ht.rsid, ht.af_UKBB, ht.ld])
    ht.write(out_prefix + ".ht")
    ht.flatten().export(out_prefix + ".txt.gz")



if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--hapmap', default=None, help='Path to hapmap file')
    parser.add_argument('--out_prefix', default=None, help='Path prefix for output dataset')
    args = parser.parse_args()

    main(args)

