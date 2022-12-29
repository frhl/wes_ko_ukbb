
library(MASS) # for glm.nb
library(pscl) # for zinf
library(ggplot2)
library(data.table)
library(argparse)

# holds all the functions to load knockouts
source("scripts/post_hoc/utils.R")

main <- function(args){

    # load gtex
    gtex <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gtex/GTEx.tstat.tsv")
    gtex <- cbind(ensembl_gene_id = gtex$ENSGID, apply(gtex[,-1], 2, function(x) x > quantile(x, probs = 0.9)))
    colnames(gtex) <- tolower(gsub("_\\s*\\([^\\)]+\\)","",colnames(gtex)))

    # categories to merge on
    gtex_categories <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gtex/GTEX.tstat.categories.genoppi.csv")
    gtex_categories$id <- tolower(gsub("_\\s*\\.[^\\.]+\\.","",gtex_categories$Tissue.genoppi))
    gtex_categories$id <- gsub("[0-9]\\.$","",gtex_categories$id)
    gtex_categories$id <- gsub("\\.","-",gtex_categories$id)
    stopifnot(sum(!gtex_categories$id %in% colnames(gtex)) == 0)

    # setup categories and annotations
    annotations <- c("pLoF_damaging_missense", "pLoF", "damaging_missense", "other_missense", "synonymous")
    categories <- colnames(gtex)[-1]

    # load data and set up nb
    dt <- fread(args$consolidated)
    model <- as.formula("n~x+norm_length+norm_length2")
    lst_anno <- lapply(annotations, function(anno){
        # subset to right annotaiton
        full <- dt[dt$annotation == anno,]
        lst_data <- list(homzygotes = full[full$subset == "hom+cis",], chets = full[full$subset == "chet+cis"])
        lst_subset <- lapply(names(lst_data), function(dt_name){
            # create subset to be used
            dt_fit <- lst_data[[dt_name]]
            dt_fit <- merge(dt_fit, gtex, by = c("ensembl_gene_id"), all.x = TRUE)
            # iterate through every category
            lst_category <- lapply(categories, function(category){
                print(paste(category, anno, dt_name))
                # remove rows that are NA
                na_rows <- rowSums(is.na(dt_fit[,c('n','norm_length','norm_length2')]))
                dt_fit <- dt_fit[!na_rows,]
                dt_fit$x <- dt_fit[[category]]
                # set up negative binomial
                x <- glm.nb(model, data = dt_fit, init.theta = 1.041, link = "log")
                fit <- data.frame(coef(summary(x)))
                colnames(fit) <- c("est", "error", "z", "p")
                # get confidence intervals
                conf <- suppressMessages(exp(cbind(coef(x), confint(x))))
                fit$ci_est <- conf[,1]
                fit$ci_lower <- conf[,2]
                fit$ci_upper <- conf[,3]
                fit$keep <- as.logical(c(0, 1, 0, 0))
                fit$id <- category
                fit$annotation <- anno
                fit$subset <- dt_name
                rownames(fit) <- NULL
                return(fit)
            })
            return(do.call(rbind, lst_category))
        })
        return(do.call(rbind, lst_subset))
    })
    # combine model
    final <- do.call(rbind, lst_anno)
    final <- final[final$keep==1,]
    fits <- merge(final, gtex_categories, all.x = TRUE)
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(fits, outfile)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--consolidated", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


