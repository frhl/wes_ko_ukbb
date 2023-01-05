
library(argparse)
library(data.table)

list_files_saige <- function(cond="none", prs="include", regex = "\\.txt\\.gz"){

    # set up paths
    wd <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb"
    step2_dir <- file.path(wd, "data/saige/output/binary/step2/min_mac4")
    step2_common_dir <- file.path(wd, "data/saige/output/binary/step2_common/min_mac4")
    step2_rare_dir <- file.path(wd, "data/saige/output/binary/step2_rare_cond/min_mac4")
    step2_combined_dir <- file.path(wd, "data/saige/output/binary/step2_rare_cond/min_mac4")
    
    # subset paths based on condions
    if (cond %in% "none"){
        files <- list.files(step2_dir, full.names = TRUE, pattern = regex)
    } else if (cond %in% "common"){
        files <- list.files(step2_common_dir, full.names = TRUE, pattern = regex)
    } else if (cond %in% "rare"){
        files <- list.files(step2_rare_dir, full.names = TRUE, pattern = regex)
    } else if (cond %in% "combined"){
        files <- list.files(step2_combined_dir, full.names = TRUE, pattern = regex)
    } else {
        stop(paste(cond, "is not a valid. Try 'none','common','rare' or 'combined'"))
    }
    
    
    # perform final subset by PRS
    files <- sort(files)
    is_prs <- grepl("locoprs.txt.gz", files)
    if (prs %in% "include"){
        return(files)
    } else if (prs %in% "exclude"){
        # exclude any PRS files
        return(files[!is_prs])
    } else if (prs %in% "only") {
        # only include PRS
        return(files[is_prs])
    } else if (prs %in% "prefer") {
        # if two files exists, one with PRS and without
        # this option will preferentially select PRS
        directory <- unique(dirname(files))
        stopifnot(length(directory) == 1)
        file_df <- data.frame(bname=basename(files))
        file_df$sans_ext <- gsub(".txt.gz", "", file_df$bname)
        file_df$sans_locoprs <- gsub("_locoprs", "", file_df$sans_ext)
        # count how many occours twice
        counts <- data.frame(table(file_df$sans_locoprs))
        colnames(counts) <- c("sans_locoprs", "n")
        file_df <- merge(file_df, counts)
        # subset to files with locoprs if there are matches
        file_df_n1 <- file_df[file_df$n == 1,]
        file_df_n2 <- file_df[file_df$n >= 2,]
        file_df_n2 <- file_df[grepl("_locoprs", file_df$bname), ]
        file_df <- rbind(file_df_n1, file_df_n2)
        files <- file.path(directory, file_df$bname)  
        return(files)
    } else {
        stop(paste(prs, "is not valid. Must be either 'include','exclude','prefer' or 'only'."))
    }
    
}


main <- function(args){

    args$p_cutoff <- as.numeric(args$p_cutoff)
    files <- list_files_saige(args$cond_step, prs=args$prs)
    d <- rbindlist(lapply(files, function(f){
        trait <- basename(f)
        trait <- stringr::str_extract(trait, "200k_.+pLoF_damaging_missense")
        trait <- gsub("200k_", "", trait)
        trait <- gsub("_pLoF_damaging_missense", "", trait)
        d <- fread(f)
        d <- d[grepl("ENSG", d$MarkerID),]
        if ("p.value_c" %in% colnames(d)) {
            d$p.value <- d$p.value_c 
            d$cond <- TRUE
        } else {
            d$cond <- FALSE
        }
        d <- d[d$p.value < args$p_cutoff,]
        if (nrow(d)){
            out <- data.table(
                trait = trait,
                chromosome = d$CHR,
                gene = d$MarkerID,
                pvalues = d$p.value,
                cond = d$cond
            )
            return(out)
        }
    }))
    # write files
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(d, outfile, sep = "\t", quote=FALSE)

}

parser <- ArgumentParser()
parser$add_argument("--cond_step", default=NULL, help = "")
parser$add_argument("--prs", default=NULL, help = "")
parser$add_argument("--p_cutoff", default=5e-7, help = "")
parser$add_argument("--out_prefix", default=NULL, help = "")
args <- parser$parse_args()

main(args)


