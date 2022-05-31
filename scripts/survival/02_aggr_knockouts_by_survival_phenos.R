library(argparse)
library(data.table)

# labels to be used
labels <- c(
    "Heterozygote" = "HET",
    "Homozygote" = "HOM",
    "Compound heterozygote" = "CHET_TRANS",
    "Compound heterozygote (cis)" = "CHET_CIS",
    "Possible Compound heterozygote" = "CHET_UNKNOWN"
)

# bridge for mapping esngid -> hgnc
bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensg_to_hgnc <- as.vector(bridge$hgnc_symbol)
names(ensg_to_hgnc) <- bridge$ensembl_gene_id



main <- function(args){

    # read on all knockouts
    ko <- do.call(rbind, lapply(1:22, function(chr){fread(gsub("CHR",chr,args$ko_file))}))
    ko <- ko[,c('gene_id',"s","knockout")]
    ko$knockout <- labels[ko$knockout]

    # read in longitidunal phenotypes
    phenotypes <- fread(args$phenotypes)
    discard_cols <- c("age_at_first_GP_record", "age_at_last_GP_record", 
                      "age_at_first_HES_record", "age_at_last_HES_record")
    phenotypes <- phenotypes[,!colnames(phenotypes) %in% discard_cols, with = FALSE]

    # also get original phenotypes and subset to NFE
    original_phenotypes <- fread(args$original_phenotypes)
    keep_samples <- original_phenotypes$eid[original_phenotypes$genetic.eur.no.fin.oct2021]
    phenotypes <- phenotypes[phenotypes$eid %in% keep_samples,]
    phenos <- colnames(phenotypes)[-1]

    # merge knockouts and phenotypes
    ko <- merge(ko, phenotypes, by.x = "s", by.y = "eid", all.x = TRUE)

        by_samples <- do.call(rbind, lapply(phenos, function(p){

        # count number of samples that are knockouts
        ko_tmp <- ko[,c('s','knockout',p), with = FALSE]
        ko_tmp[[p]] <- ifelse(!is.na(ko_tmp[[p]]), TRUE, FALSE)
        f <- as.formula(paste0("s~knockout"))
        ko_tmp <- dcast(f, data = ko_tmp, fun.aggregate = sum, value.var = p)
        ko_tmp[is.na(ko_tmp)] <- 0
        names <- colnames(ko_tmp)[-1]
        mat <- data.frame(t(matrix(colSums(ko_tmp[,-1]))))
        colnames(mat) <- names

        # annotate with all
        all_phenotypes <- phenotypes[[p]]
        mat$phenos_na <- sum(is.na(all_phenotypes))
        mat$cases_all <- sum(all_phenotypes == 1, na.rm = TRUE)
        mat$controls_all <- sum(all_phenotypes == 0, na.rm = TRUE)

        # clean up and move phenotype to first column
        mat$phenotype <- p
        mat <- mat[,c(ncol(mat), 1:(ncol(mat)-1))]
        return(mat)
    }))

    # write results
    outfile_samples <- paste0(args$out_prefix, "_by_phenotypes.txt.gz")
    fwrite(by_samples, outfile_samples, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--ko_file", default=NULL, required = TRUE, help = "")
parser$add_argument("--phenotypes", default=NULL, required = TRUE, help = "")
parser$add_argument("--original_phenotypes", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)











