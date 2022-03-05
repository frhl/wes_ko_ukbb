devtools::load_all("utils/modules/R/phasingtools")

main <- function(args){

  print(args)
  stopifnot(dir.exists(args$in_dir))
  stopifnot(dir.exists(args$out_prefix))
  stopifnot(args$alpha > 0) 
 
  files <- list.files(args$in_dir, pattern = args$in_ext, full.names = TRUE)
  if (!is.null(grep)) files <- files[grepl(args$grep, basename(files))]
  stopifnot(length(files) > 0)

  # simple analysis of non-hacked BCFtools
  #if (!is.null(args$simple)){
  # 
  #  M <- do.call(rbind, lapply(files, function(f){
  #      d <- fread(f) 
  #      d <- summarize_bcftools_trio_stats(d)  
  #      d$chr <- str_extract(f, 'chr[0-9]+')
  #      return(d)
  #  }))
  #
  #  return(M)
  #}

  # load genes generated with Biomart (GRCh38)
  genes <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/protein_coding_genes.tsv')
  stopifnot(file.exists(genes))
  genes <- genes[genes$chromosome_name %in% 1:22]
  
  # Site-specific analysis
  dt <- setDT(do.call(rbind, lapply(files, function(f){
      d <- fread(f)
      if ('MAF' %in% colnames(d)){
          
          # get data by site
          d <- aggr_ser_by_site(d)
          d$name <- tools::file_path_sans_ext(basename(f))
          
          # perform subset on chromosome and only get genes within the subset
          chr <- as.numeric(gsub("chr","",str_extract(f, 'chr[0-9]+')))
          genes_subset <- genes[genes$chromosome_name %in% chr,]
          d$chr <- chr
          
          # append gene information
          d$gene <- unlist(lapply(d$POS, function(p){
             out <- genes_subset[p > genes_subset$start_position & 
                                 p < genes_subset$end_position, ]
             return(paste0(out$hgnc_symbol, collapse = ';'))
          }))
          
          d$n_switch_cumsum <- cumsum(d$n_switch)
          d$n_tested_cumsum <- cumsum(d$n_tested)
          
          # get site processed file
          return(d)
      }
  })))
  
  # get summary
  counts <- cbind(
    data.frame(aggregate(n_tested ~ name, data = dt, FUN = sum)),
    data.frame(n_switch = aggregate(n_switch ~ name, data = dt, FUN = sum)$n_switch)
  )

  # get binominal confidence intervals for switch error estimate
  binom_ci <- do.call(rbind, lapply(1:nrow(counts), function(i){
    bconf <- binconf(
        counts$n_switch[i], 
        counts$n_tested[i], 
        alpha = args$alpha
    )
    est <- bconf[1]
    upper <- bconf[3]
    err <- upper - est
    data.frame(ser_est = est, ser_est_err = err)
  }))
 
  # Combine with original counts
  counts <- cbind(counts, binom_ci)
  counts$chr<- str_extract(counts$name, 'chr[0-9]+')
 
  # aggregate switch errors by gene (to determine how many are 
  # completely resolved without any phasing errors)
  datasets <- unique(dt$name)
  phased_genes <- do.call(rbind, lapply(datasets, function(name){
  
    bool <- dt$name == name
    d_subset <- dt[bool,]
    # aggregate switches across genes
    daggr <- aggregate(n_switch ~ gene, data = d_subset, FUN=sum)
    daggr$name <- name
    
    # get frequency of SERs per gene
    n_genes <- length(unique(daggr$gene))
    aggr_gene_freq <- as.data.frame(table(daggr$n_switch, daggr$name))
    colnames(aggr_gene_freq)[1:3] <- c('switch_errors','dataset','n')
    aggr_gene_freq$total_genes <- n_genes
    aggr_gene_freq$pct <- aggr_gene_freq$n / n_genes
      
    return(aggr_gene_freq)
  }))

  # write out stats for aggregated genes
  outfile = paste0(args$out_prefix, '_phased_genes.txt.gz')
  fwrite(phased_genes, outfile, sep = '\t')

  # append information by counts
  counts$genes_no_errors <- phased_genes[phased_genes$switch_errors == 0, ]$n
  counts$genes_total <- phased_genes[phased_genes$switch_errors == 0, ]$total_genes
  counts$genes_no_errors_pct <- phased_genes[phased_genes$switch_errors == 0, ]$pct
 
  # stats for counts overview
  outfile = paste0(args$out_prefix, '_overview.txt.gz')
  fwrite(maf_error, outfile, sep = '\t')

  # create MAF bins
  cuts <- c(0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
  dt$maf_bin <- cut(dt$MAF, breaks = cuts, levels = as.character(cuts))

  maf_error <- do.call(rbind, lapply(unique(dt$name), function(name){

     bool = dt$name %in% name
     d_subset <- dt[bool,]
     # aggregate my Minor Allele frequency Bim
     tmp <- cbind(
            data.frame(name = name),
            data.frame(aggregate(n_tested ~ maf_bin, data = d_subset, FUN = sum)),
            data.frame(n_switch = aggregate(n_switch ~ maf_bin, data = d_subset, FUN = sum)$n_switch)
     )

     # get binomial confidence interval
     bconf <- binconf(tmp$n_switch, tmp$n_tested, alpha = args$alpha)
     est <- bconf[1,]
     upper <- bconf[3,]
     err <- abs(upper - est)
     tmp$ser_est <- est
     tmp$ser_est_err <- err
     return(tmp)

  }))

  # switch errors by maf
  outfile = paste0(args$out_prefix, '_dataset_maf.txt.gz')
  fwrite(maf_error, outfile, sep = '\t')

  # aggregate my Minor Allele frequency Bim
  maf_bin_error <- cbind(
            data.frame(aggregate(n_tested ~ maf_bin, data = dt, FUN = sum)),
            data.frame(n_switch = aggregate(n_switch ~ maf_bin, data = dt, FUN = sum)$n_switch)
     )

  # get binomial confidence interval
  bconf <- binconf(maf_bin_error$n_switch, maf_bin_error$n_tested)
  est <- bconf[1,]
  upper <- bconf[3,]
  err <- abs(upper - est)
  maf_bin_error$ser_est <- est
  maf_bin_error$ser_est_err <- err
  
  # combined maf bin errors
  outfile = paste0(args$out_prefix, '_maf.txt.gz')
  fwrite(maf_error, outfile, sep = '\t')

  # full dataset used to generated the above stats
  outfile = paste0(args$out_prefix, '.txt.gz')
  fwrite(dt, outfile, sep = '\t')


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--in_ext", default='*.txt', help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--grep", default=NULL, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--alpha", default=0.05, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


