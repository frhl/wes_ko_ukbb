
library(MASS) # for glm.nb
library(ggplot2)
library(data.table)
library(argparse)

# holds all the functions to load knockouts
source("scripts/post_hoc/utils.R")

main <- function(args){

    consolidated <- args$consolidated
    stopifnot(file.exists(consolidated))
    dt <- fread(consolidated)

    # load GTEx categories
    gtex <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gtex/GTEx.tstat.tsv")
    gtex <- cbind(ensembl_gene_id = gtex$ENSGID, apply(gtex[,-1], 2, function(x) x > quantile(x, probs = 0.9)))
    colnames(gtex) <- tolower(gsub("_\\s*\\([^\\)]+\\)","",colnames(gtex)))
                            
    # categories to merge on
    gtex_categories <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gtex/GTEX.tstat.categories.genoppi.csv")
    gtex_categories$id <- tolower(gsub("_\\s*\\.[^\\.]+\\.","",gtex_categories$Tissue.genoppi))
    gtex_categories$id <- gsub("[0-9]\\.$","",gtex_categories$id)
    gtex_categories$id <- gsub("\\.","-",gtex_categories$id)
    stopifnot(sum(!gtex_categories$id %in% colnames(gtex)) == 0) 

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

    # get annotations 
    annotations <- c("pLoF","damaging_missense","pLoF_damaging_missense","other_missense","synonymous")
    categories <- colnames(gtex)[-1]

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
            mrg <- merge(mrg, gtex, all.x=TRUE, by = "ensembl_gene_id")
            
            
            lst <- lapply(categories, function(category){
                
                print(paste0("--",category))
                dt_fit <- mrg
                na_rows <- rowSums(is.na(dt_fit[,c('n','length', "gc_norm")]))
                dt_fit <- dt_fit[!na_rows,]
                dt_fit$x <- dt_fit[[category]]
               
                # adjust for transcript length/gc content 
                model <- as.formula(paste0("n~x+length+gc_norm"))#,get_mut_rate_label(anno)))
                x <- glm.nb(model, data = dt_fit, link = "log",  control=glm.control(maxit=100))
                fit <- data.frame(coef(summary(x)))
                colnames(fit) <- c("est", "error", "z", "p")
                conf <- suppressMessages(exp(cbind(coef(x), confint(x))))
                
                fit$ci_est <- conf[,1]
                fit$ci_lower <- conf[,2]
                fit$ci_upper <- conf[,3]

                fit$keep <- as.logical(c(0, 1, 0, 0))
                fit$id <- category
                rownames(fit) <- NULL
                return(fit)
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
        colnames(d1)[1:8] <- paste0("hom_",colnames(d1)[1:8])
        colnames(d2)[1:8] <- paste0("chet_",colnames(d2)[1:8])
        d <- merge(d1, d2, by = "id")
        d$annotation <- anno
        
        return(d)
    })

    d <- do.call(rbind, outer_list)
    d1 <- d[,c("id","hom_est", "hom_error", "hom_z", "hom_p", "hom_ci_est", "hom_ci_lower", "hom_ci_upper", "annotation")]
    d2 <- d[,c("id","chet_est", "chet_error", "chet_z", "chet_p", "chet_ci_est", "chet_ci_lower", "chet_ci_upper", "annotation")]
    colnames(d1) <- gsub("hom_", "", colnames(d1))
    colnames(d2) <- gsub("chet_", "", colnames(d2))
    d1$dt_name <- "homs"
    d2$dt_name <- "chets"
    dout <- rbind(d1, d2)

    # get gtex categories
    dout <- merge(dout, gtex_categories, all.x = TRUE)
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(dout, outfile)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--consolidated", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


