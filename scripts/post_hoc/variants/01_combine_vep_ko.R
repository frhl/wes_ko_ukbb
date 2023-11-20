
library(data.table)
library(argparse)

# process variants on each strand
clean_hail_list <- function(x, comma_sub = ', ') {
    x <- gsub('(\\[)|(\\])|(\\")','',x)
    x <- gsub('\\,',comma_sub,x)
    return(x)
}

# extract reference target vector by reference ID 
# vectors that are also present in target_id vector
extract_by_entry <- function(ref_id, target_id, target, unlist_result = TRUE){
    stopifnot(!is.null(target_id))
    stopifnot(!is.null(target))
    #stopifnot(length(target_id) == length(target))   
    return(
        lapply(strsplit(ref_id, split = ', '), function(row){
            lst <- lapply(row, function(id){
               return(target[target_id %in% id])
            })
            if (unlist_result){
                return(unlist(lst))
            } else {
                return(lst)
            }
        })
    )
}

# append
paste_lst <- function(x, sep = ', ') unlist(lapply(x, paste, collapse = sep))


main <- function(args){
    
    chr = args$chrom
    print(chr)
    write(paste0("Reading chr",chr), stdout())
    vep <- fread(paste0('data/mt/vep/worst_csq_by_gene_canonical/ukb_eur_wes_union_calls_200k_chr',chr,'.tsv.gz'))
    dt <- fread(paste0('data/knockouts/alt/ukb_eur_wes_200k_chr',chr,'_maf0to5e-2_pLoF_damaging_missense_all.tsv.gz'))
    dt$gts <- clean_hail_list(dt$gts)
    dt$varid <- clean_hail_list(dt$varid)
    dt$revel_score <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$csqs.revel_score))
    dt$cadd_phred <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$csqs.cadd_phred))
    dt$most_severe_csqs <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$csqs.most_severe_consequence))
    dt$AF <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$info.AF))
    dt$AC <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$info.AF))
    dt$AN <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$info.AN))
    dt$gene_symbol <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$csqs.gene_symbol))
    dt$exon <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$csqs.exon))
    dt$codon <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$csqs.codons))
    dt$amino_acids <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$csqs.amino_acids))
    dt$consequence_category <- paste_lst(extract_by_entry(dt$varid, vep$varid, vep$consequence_category))
    dt$chr <- chr

    outfile = paste0(args$out_prefix, "_chr", chr,".txt.gz")
    write(paste0("writing to ", outfile), stderr())
    fwrite(dt, outfile, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


