# module purge
# conda activate rpy
# Rscript 04_knockouts.R --in_dir derived/knockouts/211111 --in_pattern maf00_50 --in_csq ptv_damaging_missense_knockouts --out_prefix derived/summary

# setup paths and libs
library(argparse)
library(data.table)
library(ggplot2)
library(dplyr)

setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')
devtools::load_all('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/ukbtools')

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, help = "Where should the results be written?")
parser$add_argument("--in_phenotypes", default=NULL, help = "Phenotypes")
args <- parser$parse_args()

stopifnot(args$in_phenotypes %in% c('binary','cts'))
stopifnot(dir.exists(dirname(args$out_prefix)))

# knockout counts
knockouts <- list.files('derived/knockouts/211111/', full.names = TRUE)

# phenotypes
pheno_binary <- unlist(strsplit(readLines('data/phenotypes/UKBB_WES200k_binary_phenotypes_header.txt'), split = '\t'))
pheno_cts <- unlist(strsplit(readLines('data/phenotypes/UKBB_WES200k_cts_phenotypes_header.txt'), split = '\t'))

# saige results
saige_binary <- list.files('data/saige/output/combined/binary/step2/211111/', full.names = TRUE)
saige_binary <- saige_binary[file.info(saige_binary)$size > 1]
saige_cts <- list.files('data/saige/output/combined/cts/step2/211111/', full.names = TRUE)
saige_cts <- saige_cts[file.info(saige_cts)$size > 1]

# thresholds
mutations <- c('ptv','ptv_damaging_missense','synonymous')
mafs <- c('00_01','01_50','00_50')

# current params selected
mafs = '00_01' # sometimes eerors out with too high MAF
ribbon_p=0.95
out_prefix = args$out_prefix

if (args$in_phenotypes == 'binary'){
    phenotypes = pheno_binary
    saige_file = saige_binary
} else {
    phenotypes = pheno_cts
    saige_file = saige_cts
}


#out_prefix <- 'derived/tables/saige/211111_wes200k_saige_merge'

# for mapping to hgnc symbols
hgnc_bridge <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/hgnc/211026_hgnc_ensgid_link.csv')
colnames(hgnc_bridge)[2] <- 'gene_id'

for (maf in mafs){
   
    write(maf,stdout())

    for (mutation in mutations){
        
        write(mutation, stdout())

        # get the right knockout files
        bool_maf = grepl(maf, knockouts)
        bool_mutation = grepl(paste0(mutation,'_knockouts'), knockouts)
        knockout_files <- knockouts[bool_mutation & bool_maf]
        
        if (length(knockout_files) > 0){
                
            # read in knockout results and add meta data
            write(paste0('loading knockouts', knockout_files[1]), stdout())
            ko_dt <- setDT(do.call(rbind, lapply(knockout_files, zcat)))
            ko_dt$csqs_category <- mutation
            ko_dt$maf <- maf
                
            # aggregate knockout dt
            dt <- data.table(table(ko_dt$gene_id, ko_dt$csqs))
            dt <- data.table::dcast(dt, V1 ~ V2, value.var = 'N')
            colnames(dt)[1] <- 'gene_id'
                
            # deal with potential extra category
            if ('CH+HO' %in% colnames(dt)){
                dt$HO <- dt$HO + dt$`CH+HO`
                dt$`CH+HO` <- NULL
            }
                
            # aggregate knockout dt
            ko_count <- data.table(table(ko_dt$gene_id, ko_dt$csqs))
            ko_count <- data.table::dcast(ko_count, V1 ~ V2, value.var = 'N')
            colnames(ko_count)[1] <- 'gene_id'
                
            # clean up haplotypes
            dt <- ko_dt[ko_dt$knockout]
            dt$hap1 <- format_from_hail(dt$phase1)
            dt$hap2 <- format_from_hail(dt$phase2)
            dt <- dt[,c('gene_id','hap1','hap2')]
            dt <- dt[!duplicated(dt)]
                
            # sort haplotypes and remove duplicates
            haplotypes_sorted <- data.table(t(apply(dt[,c('hap1','hap2')],1,sort)))
            genes_haplotypes <- cbind(dt$gene_id, haplotypes_sorted)
            colnames(genes_haplotypes) <- c('gene_id','hap1','hap2')
            genes_haplotypes <- genes_haplotypes[!duplicated(genes_haplotypes)]
                
            #outfile_haplotypes = paste0(out_prefix,'_',maf,'_',mutation,'.tsv.gz')
            #print(outfile)
                
            # combine haplotypes by gene
            combined <- genes_haplotypes
            combined$knockout_alleles <- paste0(combined$hap1, '+', combined$hap2)
            combined <- aggregate(knockout_alleles ~ gene_id, data = combined, FUN = function(x) paste(x, collapse = '; '))
        }
        
        for (phenotype in phenotypes){
            
            write(phenotype, stdout())

            # Get the right saige files
            pheno_mutation = paste0(mutation, '_', phenotype)
            bool_maf = grepl(maf, saige_file)
            bool_pheno = grepl(phenotype, saige_file)
            bool_mutation = grepl(pheno_mutation, saige_file)
            saige_files <- saige_file[bool_pheno & bool_mutation & bool_maf]
            
            if (length(saige_files) > 0){
                
                # read in saige results and add meta data
                write(paste('loading saige files',saige_files[1]),stdout())
                saige_dt <-  setDT(do.call(rbind, lapply(saige_files, fread)))
                saige_dt$phenotype <- phenotype
                saige_dt$csqs_category <- mutation
                saige_dt$maf <- maf
                
                # calculate expected versus observed p-values
                n <- nrow(saige_dt)
                saige_dt <- saige_dt[order(saige_dt$p.value),]
                saige_dt$pvalue.observed <- saige_dt$p.value
                saige_dt$pvalue.expected <- seq(1, n)/(n + 1)
                saige_dt$clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n:1, shape1 = 1:n))
                saige_dt$cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n:1, shape1 = 1:n))
                saige_dt$FDR <- stats::p.adjust(saige_dt$pvalue.observed, method = 'fdr')

                # rename certain columns
                colnames(saige_dt)[colnames(saige_dt)=='SNPID'] <- 'gene_id'
                
                # write file
                #outfile_saige = paste0(out_prefix,'_',phenotype,'_',maf,'_',mutation,'.tsv.gz')
                #print(outfile)
                
            }
            
            if ((length(saige_files) > 0) & (length(knockout_files) > 0)){
            
                # merge saige output with knockout count and with alleles involved in knockout
                mrg <- merge(saige_dt, ko_count, all.x = TRUE, by = 'gene_id')
                mrg <- merge(mrg, combined, all.x = TRUE, by = 'gene_id')
                mrg <- merge(mrg, hgnc_bridge, all.x = TRUE, by = 'gene_id')
                
                outfile_saige = paste0(out_prefix,'_',phenotype,'_',maf,'_',mutation,'.tsv')
                write(paste('writing',outfile_saige),stdout())
                fwrite(mrg, file = outfile_saige, quote = FALSE, sep = '\t', row.names = FALSE)
                
            }
            
        }
        
    }
    
}

write('all loops ran successfully!', stdout())






