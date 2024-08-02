#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){
    
    out_prefix <- args$out_prefix
    input_path <- args$input_path
    qced_samples <- args$qced_samples
    samples_to_redact <- args$samples_to_redact
    invert <- args$invert

    dt <- fread(input_path)
   
    # filter to qced samples 
    if (!is.null(qced_samples)) {
        qced_samples <- fread(qced_samples, header=FALSE)$V1
        dt <- dt[dt$eid %in% qced_samples,]
        stopifnot(nrow(dt)>50000)
    }
   
    # redact samples 
    dt$original_order <- seq_len(nrow(dt))
    samples_to_redact <- fread(samples_to_redact, header=FALSE)$V1
    stopifnot(length(samples_to_redact)>0)

    # get samples to redact, get columns to redact 
    idx_redact <- which(dt$eid %in% samples_to_redact)
    idx_keep <- which(!dt$eid %in% samples_to_redact)
  
    # switch refact/keep if invert
    if (invert) {
        write("inverting sample redaction!", stderr())
        tmp <- idx_keep
        idx_keep <- idx_redact
        idx_redact <- tmp
    }

    # split into two
    dt_redact <- dt[idx_redact,]
    dt_keep <- dt[idx_keep,]
    print(head(dt_redact$eid))

    # redact away
    stopifnot(nrow(dt_redact)>1)
    stopifnot(nrow(dt_keep)>1)
    dt_redact <- as.data.frame(dt_redact)
    cols_to_redact <- colnames(dt_redact)
    cols_to_redact <- cols_to_redact[!cols_to_redact %in% "eid"]
    # redact away..
    for (col in cols_to_redact){
        dt_redact[[col]] <- NA
    }

    # combine and re-order by original "dt" 
    dt_combined <- rbind(dt_redact, dt_keep)
    dt_combined <- dt_combined[order(dt_combined$original_order),]
    dt_combined$original_order <- NULL

    outfile <- paste0(out_prefix, ".txt")
    write(paste("writing to", outfile), stdout())
    fwrite(dt_combined, outfile, sep = "\t", quote = FALSE, na="NA")
    outsamples <- paste0(out_prefix, ".samples")
    write(paste("writing to", outsamples), stdout())
    fwrite(data.table(eid=dt_keep$eid), outsamples, sep = "\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, help = "Where should the file be written")
parser$add_argument("--input_path", default=NULL, help = "Curated phenotypes input path")
parser$add_argument("--qced_samples", default=NULL, help = "path to qced samples")
parser$add_argument("--samples_to_redact", default=NULL, help = "Samples to set to NA")
parser$add_argument("--invert", action="store_true", default=FALSE)
args <- parser$parse_args()

main(args)

