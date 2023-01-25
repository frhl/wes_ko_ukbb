
library(MASS) # for glm.nb
library(ggplot2)
library(data.table)
library(argparse)
library(genoppi)

# holds all the functions to load knockouts
source("scripts/post_hoc/utils.R")

main <- function(args){

    consolidated <- args$consolidated
    stopifnot(file.exists(consolidated))
    dt <- fread(consolidated)

    # get transcript
    transcript <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/221216_enstid_ensgid_lengths.txt.gz")
    transcript <- transcript[,c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "chromosome_name", "length")]
    transcript$length <- (transcript$length-mean(transcript$length))/sd(transcript$length)
    transcript$length2 <- transcript$length ^ 2

    # get GC-content
    dgc <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/221229_ensgid_gc_content.txt.gz")
    dgc$gc <- dgc$percentage_gene_gc_content
    dgc$gc_norm <- (dgc$gc - mean(dgc$gc))/sd(dgc$gc)
    dgc$percentage_gene_gc_content <- NULL
    dgc_transcript <- merge(transcript, dgc, all.x = TRUE, by = c("ensembl_gene_id"))

    # get mutation rate
    mut <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/mutation_rates/samocha2014.txt.gz")
    geneset <- merge(dgc_transcript, mut, by = c("ensembl_gene_id","hgnc_symbol"), all.x = TRUE)

    hgnc_to_ensembl <- get_mapping_hgnc_to_ensembl()
    # set up kegg
    kegg <- msigdb_c2_table[grepl('KEGG', msigdb_c2_table$Set.name),]
    colnames(kegg) <- c("hgnc_symbol", "geneset")
    kegg$ensembl_gene_id <- hgnc_to_ensembl[kegg$hgnc_symbol]
    rownames(kegg) <- NULL

    # set up reactome
    reactome <- msigdb_c2_table[grepl("REACTOME", msigdb_c2_table$Set.name),]
    colnames(reactome) <- c("hgnc_symbol", "geneset")
    reactome$ensembl_gene_id <- hgnc_to_ensembl[reactome$hgnc_symbol]
    rownames(reactome) <- NULL

    # combine
    kegg_reactome <- rbind(kegg, reactome)
    kegg_reactome$hgnc_symbol <- NULL

    annotations <- c("pLoF","damaging_missense","pLoF_damaging_missense","other_missense","synonymous")
    categories <- unique(kegg_reactome$geneset)

    outer_list <- lapply(annotations, function(anno){
        print(anno)
        # subset main data and get homs/cis
        dt_subset <- dt[dt$annotation == anno,]
        lst_data <- list(
            chets = dt_subset[(dt_subset$is_chet | dt_subset$is_cis),],
            homzygotes = dt_subset[(dt_subset$is_hom | dt_subset$is_cis),]
        )

        # get OR/P-values for each dataset
        lst <- lapply(names(lst_data), function(dt_name) {

            print(paste('-',dt_name))
            full <- lst_data[[dt_name]]
            aggr_ko <- aggregate(pKO ~ ensembl_gene_id + ensembl_transcript_id, data = full, FUN=sum)
            aggr_full <- aggregate(pKO ~ ensembl_gene_id + ensembl_transcript_id, data = full, FUN=length)
            aggr <- merge(aggr_ko, aggr_full, by = c("ensembl_gene_id","ensembl_transcript_id"))
            colnames(aggr) <- c("ensembl_gene_id","ensembl_transcript_id", "n","total")
            mrg <- merge(aggr, geneset, by = c("ensembl_gene_id", "ensembl_transcript_id"), all.x = TRUE)

            
            lst <- lapply(categories, function(category){
                
                print(paste0("--",category))
                genes_in_category <- unique(kegg_reactome$ensembl_gene_id[kegg_reactome$geneset == category])
                dt_fit <- mrg
                dt_fit$x <- dt_fit$ensembl_gene_id %in% genes_in_category
                na_rows <- rowSums(is.na(dt_fit[,c('n','length', "gc_norm")]))

                fit <- NULL

                # at least X genes present                
                if (sum(dt_fit$x) > 2){
                #    tryCatch({
                     
                    model <- as.formula(paste0("n~x+length+gc_norm"))#,get_mut_rate_label(anno)))
                    x <- glm.nb(model, data = dt_fit, link = "log",  control=glm.control(maxit=100))
                    fit <- data.frame(coef(summary(x)))
                    colnames(fit) <- c("est", "error", "z", "p")
                    #print(fit)
                    #conf <- suppressMessages(exp(cbind(coef(x), confint(x))))

                    #fit$ci_est <- conf[,1]
                    #fit$ci_lower <- conf[,2]
                    #fit$ci_upper <- conf[,3]
                    fit$keep <- as.logical(c(0, 1, 0, 0))
                    fit$id <- category
                    rownames(fit) <- NULL
                #  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

                return(fit)
                }
            })

            out <- do.call(rbind, lst)
            out$dt_name <- dt_name
            return(out)
        })

        # combine the data into a single data.frame
        d1 <- lst[[1]]
        d1 <- d1[d1$keep == 1, ]
        d2 <- lst[[2]]
        d2 <- d2[d2$keep == 1, ]
        d <- rbind(d1, d2)
        d$annotation <- anno
        
        return(d)
    })

    dout <- do.call(rbind, outer_list)
    dout$keep <- NULL

    # get gtex categories
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(dout, outfile)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--consolidated", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


