library(argparse)
library(data.table)

format_from_hail <- function(x) gsub('(\\[)|(\\])|(")|(\\})|(\\{)','',x)

fread_ko <- function(f, samples, labels){
    
    # read file and 
    d <- fread(f)
    d <- d[,c(2,1,4,10), with = FALSE]
    d <- d[d$s %in% samples]
    d$varid <- format_from_hail(d$varid)

    # set consequence
    d$csqs <- labels[d$knockout]
    d$knockout <- NULL

    # get unique genes and samples
    in_genes <- unique(d$gene_id)
    in_samples <- unique(d$s)
    out_samples <- samples[!samples %in% in_samples]

    # generate grid of out samples and genes to
    # ensure that every sample has the gene used
    grid <- expand.grid(out_samples, in_genes)
    colnames(grid) <- c("s","gene_id")

    # combine the data without duplicating a sample-gene pair
    new_grid <- rbind(d[,c(1,2)], grid[,c(1,2)])
    new_grid <- new_grid[!duplicated(new_grid),]
    stopifnot(nrow(new_grid)-nrow(grid)==nrow(d))
    new_grid$varid <- NA
    new_grid$csqs <- NA 

    # output data
    d <- rbind(d, new_grid)
    colnames(d)[2:3] <- c("ensembl_gene_id", "alleles")
    return(d)
}



main <- function(args){

  stopifnot(dir.exists(args$in_dir))

  # load auxillary files
  bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
  constraints <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv")

  # subset and merge
  constraints$loeuf <- constraints$oe_lof_upper 
  constraints <- merge(constraints, bridge, all.x = TRUE, by.x = "transcript", by.y = "ensembl_transcript_id")

  # create deciles
  probs <- seq(0,1, by = 0.1)
  probs_labels <- paste0(0:9, '-', 1:10)
  deciles <- quantile(na.omit(constraints$loeuf),probs = probs)
  constraints$loeuf_decile <- cut(constraints$loeuf, deciles, labels = probs_labels)

  # subset to releveant columns
  constraints <- constraints[,c("chromosome","hgnc_symbol","ensembl_gene_id","transcript","canonical","loeuf","loeuf_decile")]
  constraints <- constraints[constraints$canonical == TRUE]
  constraints <- constraints[!duplicated(constraints),]
  constraints <- constraints[!is.na(constraints$ensembl_gene_id),]
  constraints <- constraints[order(constraints$chromosome),]

  # mapping to hgnc_symbol
  ensg_to_hgnc <- as.vector(bridge$hgnc_symbol)
  names(ensg_to_hgnc) <- bridge$ensembl_gene_id
  # mapping to loeuf scores
  ensg_to_loeuf <- as.vector(constraints$loeuf)
  names(ensg_to_loeuf) <- constraints$ensembl_gene_id
  # mapping to loeuf deciles
  ensg_to_loeuf_decile <- as.vector(constraints$loeuf_decile)
  names(ensg_to_loeuf_decile) <- constraints$ensembl_gene_id

  # labels 
  labels <- c(
    "Heterozygote" = "Heterozygote",
    "Homozygote" = "Homozygote",
    "Compound heterozygote" = "CH (Trans)",
    "Compound heterozygote (cis)" = "CH (Cis)",
    "Possible Compound heterozygote" = "CH (Unknown)"
  )
    
  # get samples we want to keep
  samples <- fread("data/phenotypes/filtered_phenotypes_cts.tsv.gz")
  samples <- samples$eid[samples$genetic.eur.no.fin.oct2021]

  # get files
  fs <- list.files('data/knockouts/alt/', pattern = "pLoF_damaging_missense_all.tsv.gz", full.names = TRUE)
  stopifnot(length(fs) > 0)

  for (f in fs){
    d <- fread_ko(f, samples, labels)
    chrom <- stringr::str_extract(f, "chr[0-9]+")
    genes <- unique(d$ensensembl_gene_id)
    
    # table containing all mutations for affected indiviudals
    d_minimal <- d[!is.na(d$csqs),]
    d_minimal$hgnc_symbol <- ensg_to_hgnc[d_minimal$ensembl_gene_id]
    d_minimal$loeuf <- ensg_to_loeuf[d_minimal$ensembl_gene_id]
    d_minimal$loeuf_decile <- ensg_to_loeuf_decile[d_minimal$ensembl_gene_id]
    d_minimal <- d_minimal[,c(1,5,2,3,4,6,7)]

    # write to main folder
    outfile <- paste0(args$out_dir, "_knockouts_",chrom,".txt.gz")
    fwrite(d_minimal, outfile, sep = "\t", quote = FALSE)

    # create folder for indiviudal genes
    thedir <- paste0(args$out_dir,"/",chrom)
    dir.create(thedir)
    for (g in genes){
        d_gene <- d[d$gene_id %in% g,]
        d_gene$hgnc_symbol <-  ensg_to_hgnc[d_gene$ensembl_gene_id]
        d_gene$loeuf_decile <-  ensg_to_loeuf_decile[d_gene$ensembl_gene_id]
        if (sum(d_gene$pKO) >= 4){
            outfile = paste0(thedir,"/",g,".txt.gz")
            fwrite(d_gene, outfile, sep = '\t', quote = FALSE)
        }
    }
  } 

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Directory in which to search for knockouts")
parser$add_argument("--in_pattern", default="pLoF_damaging_missense_all.tsv.gz", required = TRUE, help = "Pattern for grepping files.")
parser$add_argument("--out_dir", default=NULL, required = TRUE, help = "Directory for which to output results (will create a file for each gene)")
args <- parser$parse_args()

main(args)











