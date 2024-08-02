#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
devtools::load_all("utils/modules/R/prstools")
source("scripts/post_hoc/utils.R")
library(argparse)
library(data.table)

main <- function(args){

    ## get counts for actual observed CHs and HOMs
    d_plof <- read_ukb_wes_kos("pLoF")

    # count Homs
    hom_counts <- data.table(table(d_plof$gene_id[d_plof$knockout %in% "Homozygote"]))
    colnames(hom_counts) <- c("gene_id", "obs_homs")
    setkeyv(hom_counts, "gene_id")

    # count CHs
    ch_counts <- data.table(table(d_plof$gene_id[d_plof$knockout %in% "Compound heterozygote"]))
    colnames(ch_counts) <- c("gene_id", "obs_ch")
    setkeyv(ch_counts, "gene_id")

    # merge CH and homs
    obs_counts <- merge(hom_counts, ch_counts, all=TRUE)
    obs_counts[is.na(obs_counts)] <- 0
    obs_counts <- obs_counts[rev(order(obs_counts$obs_homs, obs_counts$obs_ch)),]
    obs_counts$obs_homs_per_ch <- obs_counts$obs_homs / obs_counts$obs_ch

    # read pLoF AFs 
    cols_to_keep <- c("varid", "hgnc_symbol", "gene_id","AC_A1", "AC_A2", "MAC", "MAF", "consequence_category")
    plof_af <- rbindlist(lapply(1:22, function(chr){
        dt_AC <- fread(paste0("data/mt/prefilter/pp90/ukb_wes_union_calls_200k_chr",chr,".loftee.worst_csq_by_gene_canonical.pp90.maf0_005.csqs.txt.gz"))
        dt_AC <- dt_AC[dt_AC$consequence_category == "pLoF",]
        dt_AC[, AC_A1 := info.AC]
        dt_AC[, AC_A2 := info.AN-info.AC]
        dt_AC[, MAC := pmin(AC_A1, AC_A2)]
        dt_AC[, MAF := MAC / info.AN] 
        dt_AC[, gene_id := worst_csq_by_gene_canonical.gene_id]
        dt_AC[, hgnc_symbol := worst_csq_by_gene_canonical.gene_symbol]
        dt_AC <- dt_AC[,..cols_to_keep]
        return(dt_AC) 
    }))

    genes <- unique(d_plof$gene_id[d_plof$knockout %in% c("Homozygote", "Compound heterozygote")])
    print(length(genes))

    # iterate
    seeds <- as.numeric(args$seed) 
    n_samples <- as.numeric(args$n_samples)
    combined <- do.call(rbind, lapply(seeds, function(seed_id){
        
        write(paste0("Simulation with seed=",seed_id), stderr())
        set.seed(seed_id)
        
        sim <- rbindlist(lapply(genes, function(gene_id){
            
            #write(gene_id, stderr())
            
            MAF <- plof_af$MAF[plof_af$gene_id %in% gene_id] 
            rbinom_tmp <- function(MAF) {
              rbinom(n_samples, 1, MAF)
            }

            X_1 <- sapply(MAF, rbinom_tmp)
            X_2 <- sapply(MAF, rbinom_tmp)

            homs <- which(apply((X_1 + X_2) == 2, 1, any))
            chets <- ((rowSums(X_1) >= 1) & (rowSums(X_2) >= 1))
            chets <- length(setdiff(which(chets), homs))
            homs <- length(homs)
            data.table(gene_id, chets, homs, seed_id) 
        }))
        
        return(sim)
    }))

    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(combined, outfile, sep="\t", quote=FALSE)

    # aggregate homozygotes and compound hets
    #aggr_homs <- setDT(aggregate(homs~gene_id, FUN=mean, data=combined))
    #aggr_chets <- setDT(aggregate(chets~gene_id, FUN=mean, data=combined))
    #setkeyv(aggr_homs, "gene_id")
    #setkeyv(aggr_chets, "gene_id")

    # combine and summarize
    #sim_counts <- merge(aggr_homs, aggr_chets)
    #sim_counts$homs_per_chet <- sim_counts$homs / sim_counts$chets
    #sim_counts <- sim_counts[(sim_counts$homs>0) | (sim_counts$chets>0),]
    #sim_counts <- sim_counts[order(sim_counts$homs_per_chet),]
    #setkeyv(sim_counts, "gene_id")

    # combine with observed counts
    #mrg <- merge(sim_counts, obs_counts)
    #outfile <- paste0(args$out_prefix, ".txt.gz")
    #fwrite(mrg, outfile, sep="\n", quote=FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--n_samples", default=200000, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
parser$add_argument("--seed", default=NULL, help = "?")
args <- parser$parse_args()

main(args)



