#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
devtools::load_all("utils/modules/R/prstools")

library(argparse)
library(data.table)

# get simulated paramters used
get_sim_params <- function(f){
    return(
        list(
            b = gsub("(_b_)|(_pi)","",stringr::str_extract(basename(f), "_b.+pi")),
            pi = gsub("(pi_)|(_K)","",stringr::str_extract(basename(f), "pi.+K")),
            h2 = gsub("(h2_)|(_b)","",stringr::str_extract(basename(f), "h2.+b")),
            k = gsub("(K)|(_seed)","",stringr::str_extract(basename(f), "K.+seed")),
            seed = gsub("(seed)|(_chr)","",stringr::str_extract(basename(f), "seed.+chr")),
            chr = stringr::str_extract(f, "chr[0-9]+"),
            case = gsub("case_", "",stringr::str_extract(f, "case_[0-9]+"))
        )
    )
}

# get path to phenotype simulations based on the saige path (step2)
get_sim_df_path <- function(f, dir="data/simulation/phenotypes_new", 
                            prefix="ukb_wes_union_calls", 
                            type="entries", suffix=".tsv.gz"){
    stopifnot(dir.exists(dir))
    params <- get_sim_params(f)
    fname <- paste0(
        prefix, 
        "_h2_", params$h2,
        "_b_", params$b,
        "_pi_", params$pi,
        "_K", params$k,
        "_seed", params$seed,
        "_", params$chr,
        "_", params$case,
        "_", type,
        suffix
    )
    fpath <- file.path(dir, fname)
    stopifnot(file.exists(fpath))
    return(fpath)
}

# get expected/observed P-values
get_qq_df <- function(f, ribbon_p = 0.95, AC_Allele2_cutoff=NULL){
    stopifnot(file.exists(f))
    params <- get_sim_params(f)
    d <- fread(f)
    if (!is.null(AC_Allele2_cutoff)) d <- d[d$AC_Allele2 >= AC_Allele2_cutoff,]
    if (nrow(d) > 0){
        d <- d[order(d$p.value),]
        d$p.value.expt <- get_expected_p(d$p.value, na.rm = TRUE)
        n <- length(d$p.value)
        dt <- data.table(
            ensembl_gene_id = d$MarkerID, # [order(d$p.value)],
            pvalue = d$p.value,
            pvalue.observed = -log10(d$p.value),
            pvalue.expected = -log10(d$p.value.expt),
            clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n:1, shape1 = 1:n)),
            cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n:1, shape1 = 1:n)),
            # simulation parameters
            b = as.numeric(params$b),
            pi = as.numeric(params$pi),
            h2 = as.numeric(params$h2),
            k = as.numeric(params$k),
            seed = as.numeric(params$seed),
            id = as.factor(params$case),
            # saige parameters
            AC_allele2 = as.numeric(d$AC_Allele2),
            N_case_hom = as.numeric(d$N_case_hom),
            N_ctrl_hom = as.numeric(d$N_ctrl_hom),
            fname = basename(f)
        )
        return(dt)
    }
}



main <- function(args){

    print(args)
    stopifnot(dir.exists(args$input_dir))
    stopifnot(!is.null(args$pattern))
    AC_Allele2_cutoff <- as.numeric(args$ac_allele2_cutoff)
    stopifnot(AC_Allele2_cutoff >= 0)

    # get files to run
    files <- list.files(args$input_dir, pattern=args$input_pattern)
    files <- files[!grepl("index", files)]
    print(head(files))

    # go over each file
    lst <- list()
    for (f in files){
        
        write(f, stderr())
        # get saige file and simulated pheno file
        qq_df <- get_qq_df(f)
        sim_path <- get_sim_df_path(f, 0.95, AC_Allele2_cutoff)
        sim_cmd <- paste("zcat", sim_path, "| grep -v NA")
        sim_df <- suppressMessages(fread(sim_cmd))
        colnames(sim_df)[colnames(sim_df) == "rsid"] <- "ensembl_gene_id"
        
        # get markers that an effect
        markers_with_effect <- sim_df[,c('ensembl_gene_id','theta')]
        markers_with_effect <- markers_with_effect[markers_with_effect$theta > 0]
        markers_with_effect <- markers_with_effect[!duplicated(markers_with_effect),]
        tmp_df <- merge(qq_df, markers_with_effect, by = "ensembl_gene_id", all.x = TRUE)
        tmp_df$theta[is.na(tmp_df$theta)] <- 0

        # get cases per marker
        markers_with_cases <- sim_df[,c('ensembl_gene_id','case')]
        markers_with_cases <- markers_with_cases[markers_with_cases$case == TRUE,]
        markers_with_cases <- aggregate(case~ensembl_gene_id, data=markers_with_cases, FUN=sum)
        tmp_df <- merge(tmp_df, markers_with_cases, by = "ensembl_gene_id", all.x = TRUE)
        tmp_df$case[is.na(tmp_df$case)] <- 0
        
        lst[[f]] <- tmp_df
    }

    # merge and write
    final_df <- rbindlist(lst)
    out <- paste0(args$out_prefix, ".txt.gz")
    fwrite(final_df, out, sep = "\t")

    # also run ROC
    h2s <- unique(final_df$h2)
    bs <- unique(final_df$b)
    lst <- list()
    for (cur_h2 in h2s){
        for (cur_b in bs){
            cur_df <- final_df[(final_df$h2 == cur_h2) & (final_df$b == cur_b)]
            plims <- c(min(cur_df$pvalue), max(cur_df$pvalue))
            pseq <- seq_log(from=plims[1], to=plims[2], length.out=5000)
            dt <- rbindlist(lapply(pseq, function(p){
                cur_df_by_p <- cur_df[,]
                total <- nrow(cur_df)
                bool_true_positive <- ((cur_df_by_p$theta != 0) & (cur_df$pvalue <= p))
                bool_true_negative <- ((cur_df_by_p$theta == 0) & (cur_df$pvalue > p))
                bool_false_positive <- ((cur_df_by_p$theta == 0) & (cur_df$pvalue <= p))
                bool_false_negative <- ((cur_df_by_p$theta != 0) & (cur_df$pvalue > p))
                sensitivity <- sum(bool_true_positive) / (sum(bool_true_positive)+sum(bool_false_negative))
                specificity <- sum(bool_true_negative)/ (sum(bool_true_negative)+sum(bool_false_positive))
                return(data.table(p, sensitivity, specificity))
            }))
            dt$b <- cur_b
            dt$h2 <- cur_h2
            outname <- paste0("h2=", cur_h2,", b=",cur_b)
            dt$id <- outname
            lst[[outname]] <- dt
        }
    }
    # combine ROC data
    roc_data <- rbindlist(lst)
    roc_data$h2_label <- paste0("h2=",roc_data$h2)
    out_roc <- paste0(args$out_prefix, ".roc.txt.gz")
    fwrite(roc_data, out_roc, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_dir", default=NULL, help = "?")
parser$add_argument("--ac_allele2_cutoff", default=0, help = "?")
parser$add_argument("--input_pattern", default=NULL, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

