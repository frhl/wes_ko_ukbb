
source("scripts/post_hoc/utils.R")
library(data.table)
library(argparse)
library(stringr)

main <- function(args){

    chrom = args$chrom
    annotation = args$annotation
    out_prefix = args$out_prefix
    path_phenotypes = args$path_phenotypes
    path_header = args$path_header

    # iterate over phenotypes
    df_phenotype <- fread(path_phenotypes)
    phenotypes <- readLines(path_header)

    # read knockouts
    d <- read_ukb_wes_kos(annotation = annotation, chromosomes = chrom)
    types_to_keep <- c("Compound heterozygote (cis)", "Compound heterozygote", "Homozygote")
    d <- d[d$knockout %in% types_to_keep, ]
    d$is_cis <- d$knockout %in% "Compound heterozygote (cis)"
    d$is_chet <- d$knockout %in% "Compound heterozygote"
    d$is_hom <- d$knockout %in% "Homozygote"
    d <- d[(d$is_cis) | (d$is_chet) | (d$is_hom), ]
    genes <- unique(d$gene_id)
    write(paste("Running", length(genes), "genes.."), stderr())

    # search for varaints in our format e.g. chr21:41871518:G:A;chr21:41878840:C:CG
    regexify_varid <- function(v) return(paste0("(^",v,"$)|(^",v,";)|(;",v,";)|(;",v,"$)"))
    # remove duplicate unordered rows in matrix format
    duplicated_varids <- function(x) duplicated(do.call(rbind, lapply(1:nrow(x), function(idx) sort(data.frame(t(x[idx,]))[,1]) )))

    genes <- unique(d$gene_id)

    out <- do.call(rbind, lapply(phenotypes, function(pheno){   
    #for (pheno in phenotypes) {
        case_eid <- df_phenotype$eid[df_phenotype[[pheno]]]
        d_pheno <- d[d$s %in% case_eid,]
        write(paste("-", pheno, chrom), stderr())
        if (nrow(d_pheno) > 0) {
            out <- do.call(rbind, lapply(genes, function(g){
                d_gene <- d_pheno[d_pheno$gene_id %in% g,]
                if (nrow(d_gene) > 0){
                    write(paste("--", g, chrom), stderr())
                    haplotypes <- c("chet","hom","cis") #, "opposite")
                    out_hap <- do.call(rbind, lapply(haplotypes, function(h){
                        if (h == "chet"){
                            d_gene_hap <- d_gene[d_gene$is_chet | d_gene$is_hom,]
                        } else if (h == "hom"){
                            d_gene_hap <- d_gene[d_gene$is_hom,]
                        } else if (h == "cis"){
                            d_gene_hap <- d_gene[d_gene$is_cis,]
                        }
                        # iterate over variants either co-occuring on the same
                        # haplotype or oppoposite haplotypes
                        if (nrow(d_gene_hap) > 0){
                            # dont want to visit the same variant twice
                            varids <- unique(unlist(strsplit(d_gene_hap$varid, split = ";")))
                            grid <- expand.grid(varids, varids)
                            grid <- grid[!duplicated_varids(grid),]
                            colnames(grid) <- c("v1", "v2")
                            do.call(rbind, lapply(grid$v1, function(v1){
                                re_v1 <- regexify_varid(v1)
                                v1_in_d <- grepl(re_v1, d_gene_hap$varid)
                                do.call(rbind, lapply(grid$v2, function(v2){    
                                    # if the variant is the same it must
                                    # be the same row being counted or a homozygote
                                    if ((v1 != v2) | (h == "hom")){
                                        re_v2 <- regexify_varid(v2)
                                        v2_in_d <- grepl(re_v2, d_gene_hap$varid)
                                        occ <- data.table(table(v2_in_d, v1_in_d))
                                        occ$v1 <- v1 
                                        occ$v2 <- v2
                                        occ$g <- g
                                        occ$chrom <- chrom
                                        occ$haplotype <- h
                                        occ$phenotype <- pheno
                                        # only keep those where both variants
                                        # are present
                                        occ <- occ[(occ$v1_in_d == TRUE) &
                                                   (occ$v2_in_d == TRUE)]
                                        occ$v2_in_d <- NULL
                                        occ$v1_in_d <- NULL
                                        return(occ)
                                    }
                                }))
                            }))
                        }
                 
                    }))
               return(out_hap)
               } 
           }))
       }
   }))
   out <- out[out$N > 0,]
   outfile <- paste0(out_prefix,".txt.gz")
   write(paste("writing", outfile), stderr()) 
   fwrite(out, outfile, sep = "\t")
     

}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "")
parser$add_argument("--annotation", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_phenotypes", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_header", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


