# module purge
# conda activate rpy

# setup paths and libs
library(argparse)
library(data.table)
library(ggplot2)
library(dplyr)

setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')
devtools::load_all('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/ukbtools')

main <- function(args){

  # get geneomic coordinates and max contig length
  contig_lens <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/hail/contigs/GRCh38_lengths.txt')
  contig_lens$chromosome <- gsub('chr','',contig_lens$chromosome)
  positions <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/211124_ensgid_to_grch38_pos.tsv.gz')
  colnames(positions)[1] <- 'gene_id'
  positions <- merge(positions, contig_lens, all.x = TRUE, by.x = 'chromosome_name', by.y = 'chromosome')
  print(str(positions))
  print(str(contig_lens))

  # Extract intervals
  in_dir = args$in_dir
  in_prefix = args$in_prefix
  out_prefix = args$out_prefix

  # check that files exists
  files <- list.files(in_dir, pattern = paste0('^',in_prefix), full.names = TRUE)
  stopifnot(length(files) > 0)

  # thresholds
  stopifnot(file.exists(args$in_phenotypes))
  phenotypes <- unlist(strsplit(readLines(args$in_phenotypes), split = '\t'))
  #phenotypes <- unlist(strsplit(readLines('data/phenotypes/UKBB_WES200k_binary_phenotypes_header.txt'), split = '\t'))
  mutations <- c('ptv','ptv_damaging_missense','synonymous')
  mafs <- c('00_01')
  padding=0

  gene_list <- list()

  for (maf in mafs){
      
      for (mutation in mutations){
          
          for (phenotype in phenotypes){
              
              # get the right file
              id <- paste0('merge_',phenotype,'_',maf,'_',mutation,'.tsv')
              bool_id <- grepl(id, files)
              if (sum(bool_id) > 0){
                 
                  # read in the file and subset by FDR 
                  x <- files[bool_id]
                  dt <- fread(x)
                  dt <- dt[dt$FDR <= 0.1,]

                  # merge with positions and write out 
                  if (nrow(dt) > 0){
                      mrg <- merge(dt, positions, by = 'gene_id',all.x = TRUE)
                      out <- data.table(
                          gene = mrg$gene_id,
                          FDR = mrg$FDR,
                          maf = maf,
                          mutation = mutation,
                          phenotype = phenotype, 
                          contig = mrg$chromosome_name,
                          start = mrg$start_position,
                          end = mrg$end_position,
                          contig_length = mrg$length
                      )
                      gene_list[[id]] <- out

                  }
                  
              }
   
          }
          
      }
      
  }


  # write result
  d_out <- do.call(rbind, gene_list)
  outfile <- paste0(out_prefix, '_saige_sig_genes_intervals.txt')
  fwrite(d_out, file = outfile, row.names = FALSE, sep = '\t')
  write('intervals were written successfully',stdout())

}

parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, help = "Where should the results be written?")
parser$add_argument("--in_dir", default=NULL, help = "in directory")
parser$add_argument("--in_prefix", default=NULL, help = "prefix of files to be loaded")
parser$add_argument("--in_phenotypes", default=NULL, help = "Path to list of phenotypes to be loaded.")
args <- parser$parse_args()


