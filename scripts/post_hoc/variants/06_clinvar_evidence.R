
source("scripts/post_hoc/utils.R")
library(data.table)
library(argparse)
library(stringr)

main <- function(args){
   
    # get knockouts
    kos <- read_ukb_wes_kos(annotation = "pLoF_damaging_missense")
    colnames(kos)[colnames(kos) == "gene_id"] <- "ensembl_gene_id"
    kos <- kos[kos$pKO >= 0.5,]

    # only extract chet/hom status if the phenotype 
    # is defined (i.e. not NA) for a given sample.
    phenotype_df <- fread(args$path_phenotypes)
    d <- fread(args$path_hits_to_analyze)

    # for now, we need to do slightly different pre-processing
    if (args$path_hits_type %in% "saige") {
        # load hits that we want to get clinvar data on
        d$phenotype <- gsub("chr[0-9]+\\_","",d$phenotype)
        d <- d[d$p.value < (0.05 / (313 * 1143)),]
        d <- d[d$N_ko_case >= 2]
        d <- d[,1:3]
        d <- d[!duplicated(d),]    
    } else if (args$path_hits_type %in% "hazard") {
        # so far no post-processing needed 
        d <- d
    }

    # ensure right columsn
    stopifnot("MarkerID" %in% colnames(d))
    stopifnot("CHR" %in% colnames(d))
    stopifnot("phenotype" %in% colnames(d))
    #stopifnot(any(d$phenotype %in% colnames(phenotype_df)))

    # go over all hits and analyze
    for (idx in 1:nrow(d)){
        gene <- d$MarkerID[idx]
        chr <- d$CHR[idx]
        phenotype <- d$phenotype[idx]
        eids <- phenotype_df$eid[phenotype_df[[phenotype]]]
        path_info_chr <- gsub("chrCHR", chr, args$path_info)
        path_clinvar_chr <- gsub("chrCHR", chr, args$path_clinvar)
        write(paste("running", gene, phenotype), stderr())
        evidence <- summarize_variant_evidence(
            gene = gene,
            phenotype = phenotype,
            path_info = path_info_chr,
            path_clinvar = path_clinvar_chr,
            kos = kos,
            eids = eids
        )
        outfile <- paste0(args$out_prefix,"_",phenotype,"_",gene,".txt")
        fwrite(evidence, outfile, sep = "\t")
    } 


}

summarize_variant_evidence <- function(gene, phenotype, path_info, path_clinvar, kos, eids = NULL){
    
    stopifnot(nrow(kos)>0)
    stopifnot(nrow(phenotype) > 0)
    stopifnot(file.exists(path_info))
    stopifnot(file.exists(path_clinvar))
    
    # check same chromosome has been used
    chr_info <- stringr::str_extract(path_info, "chr[0-9]+")
    chr_clinvar <- stringr::str_extract(path_clinvar, "chr[0-9]+")
    stopifnot(chr_info == chr_clinvar)
    
    # load sames that are 
    if (is.null(eids)) {
        eids <- unique(kos$s) 
    }
    stopifnot(length(eids) > 1)
    samples <- kos
    samples_with_genes <- samples[samples$ensembl_gene_id %in% gene,]
    variants <- unique(unlist(strsplit(samples_with_genes$varid, split = ";")))
    loci <- stringr::str_extract(variants, "chr[0-9]+\\:[0-9]+")
    
    # only keep variants with the consequences
    csqs_to_keep <- c("damaging_missense","pLoF")

    # get clinvar data
    cv <- fread(path_clinvar)
    colnames(cv) <- gsub("worst_csq_by_gene_canonical\\.","",colnames(cv))

    # exported as hail locus/alleles and need to reformat
    cv <- cv[cv$locus %in% loci,]
    cv$alleles <- gsub('(\\")|(\\[)|(\\])', '', cv$alleles)
    cv$alleles <- gsub(",","\\:",cv$alleles)
    cv$varid <- paste0(cv$locus, ":", cv$alleles)
    stopifnot(nrow(cv)>0) # standard import
    
    # how many times does a variant participate in hom/chet ---------
    lst <- as.list(variants)
    lst <- lapply(lst, function(x) list(chet=0, hom=0))
    names(lst) <- variants
    # get how many times a variants participates in a chet/hom
    for (idx in 1:nrow(samples_with_genes)){
        vs <- unlist(strsplit(samples_with_genes$varid[idx], split = ";"))
        for (v in vs){
            if (length(vs) > 1) {
                lst[[v]][["chet"]] <- lst[[v]][["chet"]] + 1
            } else {
                lst[[v]][["hom"]] <- lst[[v]][["hom"]] + 1
            }
        }
    }
    # get how many times a variant participates in a hom/chet
    lst_combined <- rbindlist(lst)
    colnames(lst_combined) <- paste0("AC_",colnames(lst_combined))
    lst_combined$varid <- names(lst)
                  
    # merge with participates in
    if (! all(cv$varid %in% lst_combined$varid)) warning(paste(gene,":",phenotype, "some variants not in clinvar."))
    stopifnot(all(lst_combined$varid %in% cv$varid))
    cv <- merge(cv, lst_combined, by="varid", all=TRUE)
    stopifnot(nrow(cv) > 0 ) # chet/hom
                  
    # get actual AC by gene -------
    info <- fread(path_info)
    info <- info[info$consequence_category %in% csqs_to_keep,]
    info <- info[info$varid %in% cv$varid,]
    cols_to_keep <- c("varid","info.AF","info.AC")
    info <- info[,..cols_to_keep]
    colnames(info) <- gsub("info\\.","",colnames(info))

    # merge with participates in
    if (! all(cv$varid %in% info$varid)) warning(paste(gene,":",phenotype, "some clinvar variants not in info."))
    stopifnot(all(info$varid %in% cv$varid))
    cv <- merge(cv, info, by="varid", all=TRUE)
    stopifnot(nrow(cv) > 0 ) # info
                  
    # process a few columns to make them look nicer
    cv$clinvar.CLNDN <- gsub("\\|",";",gsub('(\\")|(\\[)|(\\])', '', cv$clinvar.CLNDN))
    cv$clinvar.CLNHGVS <- gsub("\\|",";",gsub('(\\")|(\\[)|(\\])', '', cv$clinvar.CLNHGVS))
    cv$clinvar.CLNREVSTAT <- gsub(",",";",gsub('(\\")|(\\[)|(\\])', '', cv$clinvar.CLNREVSTAT))
    cv$clinvar.CLNSIG <- gsub('(\\")|(\\[)|(\\])', '', cv$clinvar.CLNSIG)
    cv$clinvar.CLNSIGCONF <- gsub("\\|",";",gsub('(\\")|(\\[)|(\\])', '', cv$clinvar.CLNSIGCONF))
                     
    # only damaging missense and pLoF
    cv <- cv[cv$consequence_category %in% csqs_to_keep,]

    # cols to keep and order
    cols_to_keep <- c(
        "gene_id",
        "transcript_id",
        "protein_id",
        "gene_symbol",
        "varid",
        "AC",
        "AF",
        "AC_chet",
        "AC_hom",
        "most_severe_consequence",
        "consequence_category",
        "clinvar.CLNVC",
        "clinvar.CLNSIG",
        "ccds",
        "amino_acids",
        "codons",
        "exon",
        "lof",
        "lof_flags",
        "revel_score",
        "cadd_phred",
        "polyphen_score",
        "polyphen_prediction",
        "sift_prediction",
        "clinvar.CLNHGVS",
        "clinvar.CLNSIGCONF",
        "clinvar.CLNDN"
    )

    return(cv[,..cols_to_keep])
    
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_clinvar", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_info", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_phenotypes", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_hits_to_analyze", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_hits_type", default="saige", required = TRUE, help = "Either 'saige' or 'hazard'.")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


