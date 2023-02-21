#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
devtools::load_all("utils/modules/R/prstools")

library(argparse)
library(data.table)

# get corresponding simulation file
grab_sim_file <- function(file, simdir = "data/simulation/phenotypes/", hail = "cols"){
    f <- gsub("(_pLoF_damaging_missense.txt)|(case_)","",basename(file))
    f <- paste0(simdir,f, "_",hail,".tsv.gz")
    if (!file.exists(f)) stop(paste(f, "could not be found!"))
    return(f)
}

# get corresponding file with entries
get_entry_file <- function(x){
    directory <- "data/simulation/phenotypes/"
    y <- gsub("(_case)","",x)
    y <- gsub("\\.txt", "", y)
    final <- paste0(directory,y,'_entries.tsv.gz')
    return(final)
}

# process saige df and process with upper/lower CIs
get_qq_df2 <- function(files, get_sim_file = FALSE){
    ribbon_p <- 0.95
    d <- do.call(rbind, lapply(files, function(f){
        stopifnot(file.exists(f))
        d <- fread(f)
        
        extracted_pi <- gsub("(pi_)|(_K)","",stringr::str_extract(basename(f), "pi.+K"))
        extracted_vars <- gsub("(var_)|(_pi)","",stringr::str_extract(basename(f), "var.+pi"))
        pis <- unlist(strsplit(extracted_pi, split = '_'))
        vars <- unlist(strsplit(extracted_vars, split = '_'))
        
        if (nrow(d) > 0){
            
            # get qq-plot data
            d <- d[order(d$p.value),]
            # need to order this and keed cases and controls to check inflation
            d$p.value.expt <- get_expected_p(d$p.value, na.rm = TRUE)
            n <- length(d$p.value)
            dt <- data.table(
                ensembl_gene_id = d$MarkerID, # [order(d$p.value)],
                pvalue = d$p.value,
                pvalue.observed = -log10(d$p.value),
                pvalue.expected = -log10(d$p.value.expt),
                clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n:1, shape1 = 1:n)),
                cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n:1, shape1 = 1:n)),
                #gsub("(pLoF_damaging_missense.txt)|(ukb_eur_)|(K0.1_)|(_chr22)", "", basename(f)),
                #phenofile = gsub("_pLoF_damaging_missense.txt","_phenos.tsv.gz", basename)
                h2 = as.numeric(gsub("(h2_)|(_var)","",stringr::str_extract(basename(f), "h2.+var"))),
                seed = as.numeric(gsub("(seed)|(_chr)","",stringr::str_extract(basename(f), "seed.+chr"))),
                pi_1 = as.numeric(pis[1]),
                pi_2 = as.numeric(pis[2]),
                var_1 = as.numeric(vars[1]),
                var_2 = as.numeric(vars[2]),
                AC_allele2 = as.numeric(d$AC_Allele2),
                N_case_hom = as.numeric(d$N_case_hom),
                N_ctrl_hom = as.numeric(d$N_ctrl_hom)
            )
            
            # grab actual file used for simulation
            if (get_sim_file){
                print(f)
                simfile <- grab_sim_file(f)
                dsim <- fread(simfile)
                dt$sim_h2 <- round(var(dsim$y_no_noise_rescaled), 5)
                dt$sim_y <- round(var(dsim$y))
                dt$sim_cases <- sum(dsim$case)
                dt$sim_controls <- sum(!dsim$case)
                dt$var_by_rec <- var(dsim$y_no_noise_rec) / (var(dsim$y_no_noise_add) + var(dsim$y_no_noise_rec))
                dt$var_by_rec <- round(dt$var_by_rec, 4)
            }
            
            dt$fname <- basename(f)
            return(dt)
        } else {
            return(NULL)
        }
    }))
    d$label <- factor(paste0("sim",d$seed), levels = paste0("sim",1:100))
    return(d)
}


main <- function(args){

    print(args)
    stopifnot(dir.exists(args$input_dir))

    # get files to be aggregated
    pattern <- args$seed_regex
    files <- list.files(args$input_dir, pattern = pattern, full.names = TRUE)
    files <- files[!grepl("index", files)]
    files <- files[grepl(args$var_regex, files)]

    # load files and and links
    d <- get_qq_df2(files, get_sim_file = FALSE)

    # grab corresponding entries geneated with hail
    d$fentry <- get_entry_file(d$fname)
    file_df <- d[,c("fname","fentry")]
    file_df <- file_df[!duplicated(file_df)]
    file_df$h2 <- as.numeric(gsub("(h2_)|(_var)","",stringr::str_extract(basename(file_df$fname), "h2.+var")))
    file_df$seed <- as.numeric(gsub("(seed)|(_chr)","",stringr::str_extract(basename(file_df$fname), "seed.+chr")))
    file_df$id <- as.numeric(gsub("(case_)|(_pLoF)","",stringr::str_extract(basename(file_df$fname), "case.+pLoF")))

    # itereate over each file and grab auxillary associated files
    lst <- lapply(1:nrow(file_df), function(idx){
        write(idx, stderr())
        file_entry <- file_df$fentry[idx]
        file_saige <- file_df$fname[idx]
        command <- paste("zcat",file_entry,"| grep -v NA")
        entry_dt <- fread(cmd = command)
        
        # get knockout count
        ko_dt <- data.table(table(entry_dt$rsid, entry_dt$pKO != 0))
        ko_dt <- ko_dt[ko_dt$V2 == TRUE,]
        ko_dt <- ko_dt[,c("V1","N")]
        colnames(ko_dt) <- c("ensembl_gene_id","N")
        
        # get beta/theta estimates
        entry_dt <- entry_dt[,c("rsid","theta", "beta")]
        entry_dt <- entry_dt[!duplicated(entry_dt),]
        colnames(entry_dt)[1] <- "ensembl_gene_id"
        
        # get pvalue observed
        saige_dt <- d[d$fname %in% file_saige,]
        saige_dt <- saige_dt[,c("ensembl_gene_id","pvalue.observed")]
        saige_dt <- saige_dt[!duplicated(saige_dt),]
        dt <- merge(entry_dt, saige_dt)
        dt <- merge(dt, ko_dt, all.x = TRUE)
        dt$seed <- file_df$seed[idx]
        dt$h2 <- file_df$h2[idx]
        dt$analysis <- basename(file_saige)
        dt$id <- file_df$id[idx]
        return(dt)
    })

    # get combined data
    combined <- do.call(rbind, lst)
    combined$id <- factor(combined$id)
    combined$id <- as.numeric(gsub("\\.txt","",stringr::str_extract(combined$analysis, "[0-9]+.txt")))

    # write to file
    out <- paste0(args$out_prefix, ".txt.gz")
    fwrite(combined, out, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_dir", default=NULL, help = "?")
parser$add_argument("--seed_regex", default=NULL, help = "?")
parser$add_argument("--var_regex", default=NULL, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

