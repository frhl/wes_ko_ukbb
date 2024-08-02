#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(digest)

gsub_to_dosage <- function(gt){
    gt <- gsub("\\|", "\\\\", gt)
    gt <- gsub("0/1", "1", gt)
    gt <- gsub("1/0", "1", gt)
    gt <- gsub("0/0", "0", gt)
    gt <- gsub("1/1", "2", gt)
    return(gt)
}

hash_genotypes <- function(G, algo = "xxhash32"){
    require(digest)
    gt_string <- unlist(apply(G, 1, function(x) as.character(paste(x, collapse = '-'))))
    dosage_string <- gsub_to_dosage(gt_string)
    the_hash <- unlist(lapply(dosage_string, function(x) digest(x, algo=algo)))
    return(the_hash)
}

main <- function(args){

    stopifnot(file.exists(args$in_vcf))
    stopifnot(file.exists(args$pheno_file))
    stopifnot(file.exists(args$covariates))
    stopifnot(dir.exists(dirname(args$out_prefix)))
    
    # read pheno file
    main_pheno_df <- fread(args$pheno_file)
 
    # read header one time
    header_cmd <- paste0("zcat ", args$in_vcf, " | grep '#CHR' ")
    header <- fread(cmd = header_cmd, header = FALSE)
    header <- as.character(t(header)[,1])

    # zero index for chunk
    total_chunks <- as.numeric(args$total_chunks)
    chunk <- as.numeric(args$chunk) - 1
    lines_per_chunk <- as.numeric(args$lines_per_chunk)

    # indexes for each chunk
    idx_start <- (chunk * lines_per_chunk) + 1
    idx_end <- ((chunk + 1) * lines_per_chunk ) + 1
   
    msg <- paste("Reading chunk",chunk + 1,"of", total_chunks, "with indexes", idx_start, "to", idx_end)
    write(msg, stdout()) 
    
    # read currentn chunk avoiding header line
    cmd <- paste0("zcat ", args$in_vcf, " | grep -v '#' | sed '", idx_start, ",", idx_end, "!d'")
    d <- fread(cmd = cmd, header = FALSE)
    colnames(d) <- header

    # get columns with genotypes / metaid
    genotype_cols <- suppressWarnings(!is.na(as.numeric(colnames(d))))
    id_cols <- suppressWarnings(is.na(as.numeric(colnames(d))))
    
    # setup main data.table from which we will append stuff to
    id <- d[,id_cols, with = FALSE]
    out <- data.table(
          chr = id$`#CHROM`, 
          pos = id$POS, 
          id = id$ID, 
          ref = id$REF,
          alt = id$ALT
      )

    # copy root markers
    out_ac <- out
    out_hash <- out

    # Need to ensure that G is all numerics
    G <- d[,genotype_cols, with = FALSE]
    suppressWarnings(G[, names(G) := lapply(.SD, as.numeric)])

   # read in phenotypes
    phenotypes <- unlist(strsplit(args$phenotypes, split = ","))
    phenotypes <- gsub(" ", "", phenotypes)
    phenotypes <- phenotypes[phenotypes %in% colnames(main_pheno_df)]
      
    for (phenotype in phenotypes) {
       
      write(paste("opening",phenotype,"for", args$in_vcf, "at chunk", chunk + 1, "of", total_chunks), stdout())
      pheno_df <- main_pheno_df
      
      # are there samples with missing covariates?
      col_cov <- unlist(strsplit(readLines(args$covariates), split = ","))
      lst <- lapply(col_cov, function(col){row_ok <- is.na(pheno_df[[col]])})
      missing_cov <- rowSums(do.call(cbind, lst)) > 0
      pheno_df <- pheno_df[!missing_cov, ]
      #msg <- paste("Note: Removed", sum(missing_cov),"samples with missing covariates.")
      #write(msg, stdout())

      # get defined phenotypes
      defined_phenos <- !is.na(pheno_df[[phenotype]])
      eid_with_defined_phenos <- pheno_df$eid[defined_phenos]
      stopifnot(length(eid_with_defined_phenos) > 0)

      # get subset of G which contains defiend phenotypes
      G_with_defined_phenos <- colnames(G) %in% eid_with_defined_phenos
      G_subset <- G[,G_with_defined_phenos, with = FALSE]

      # Get allele count and hash for genotypes
      G_subset_AC <- rowSums(G_subset, na.rm = TRUE)
      G_subset_hash <- hash_genotypes(G_subset, "xxhash32")

      # generate allele count outfile
      out_ac[[phenotype]] <- G_subset_AC
      #lst_ac[[phenotype]] <- G_subset_AC
      #fwrite(out, outfile_ac, sep = '\t')

      # generate allele count outfile
      #outfile_hash <- paste0(args$out_prefix,"_",phenotype,"_hash.txt.gz")
      out_hash[[phenotype]] <- G_subset_hash
      #lst_hash[[phenotype]] <- G_subset_hash

    }

    # combine files
    out_ac <- data.table(do.call(cbind, out_ac))
    out_hash <- data.table(do.call(cbind, out_hash))

    # save to file 
    outfile_ac <- paste0(args$out_prefix,".AC.txt.gz")
    outfile_hash <- paste0(args$out_prefix,".hash.txt.gz")
    fwrite(out_ac, outfile_ac, sep = '\t')
    fwrite(out_hash, outfile_hash, sep = '\t')



}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_vcf", default=NULL, required = TRUE, help = "Path to vcf")
parser$add_argument("--vcf_lines", default=NULL, required = TRUE, help = "how lines in vcf")
parser$add_argument("--lines_per_chunk", default=1000, help = "how many lines per chunk")
parser$add_argument("--chunk", default=NULL, required = TRUE, help = "current chunk (parallel)")
parser$add_argument("--total_chunks", default=NULL, required = TRUE, help = "total number of chunks (parallel)")
parser$add_argument("--pheno_file", default=NULL, required = TRUE, help = "path to phenotype file")
parser$add_argument("--phenotypes", default=NULL, required = TRUE, help = "string of comma seperated phenotypes")
parser$add_argument("--covariates", default=NULL, required = TRUE, help = "string of comma sepearted covariates")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Out prefix path")
args <- parser$parse_args()

main(args)









