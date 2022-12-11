
library(data.table)
library(argparse)


main <- function(args){
   
    # load unrelated samples
    unrelated <- fread(args$unrelated_path)
    ph <- fread(args$phenotypes_path)
    header <- fread(args$phenotypes_header, header = FALSE)
    covars <- readLines(args$covars_path)
    covars <- unlist(strsplit(covars, split = ',')) 
    
    # subset to unrelated
    ph <- ph[ph$eid %in% unrelated$s, ]
    ph$eid <- as.character(ph$eid)    

    # load all knockouts
    knockouts_paths <- list.files(args$knockout_dir, pattern = args$knockout_pattern, full.names = TRUE)    

    # load knockouts
    lst_knockouts <- lapply(knockouts_paths, function(path){
        d <- fread(path)
        d$s <- as.character(d$s)
        d <- d[,c("gene_id","s","knockout","pKO")]
        return(d)
    })

    # names to keep
    cis <- "Compound heterozygote (cis)"
    chet <- "Compound heterozygote"
    homs <- "Homozygote"
    het <- "Heterozygote"

    # check for common knockouts
    dt_all <- do.call(rbind, lst_knockouts)
    counts <- data.table(table(dt_all$gene_id[dt_all$pKO > 0]))
    counts <- counts[rev(order(counts$N))]
    common_plofs <- counts[counts$N > 10000,]
    dt_all <- dt_all[!(dt_all$gene_id %in% common_plofs$V1),]
    print(paste("excluded",nrow(common_plofs),"common knockouts."))

    dt_het <- dt_all[dt_all$knockout %in% het,]
    dt_cis <- dt_all[dt_all$knockout %in% cis,]
    dt_ko <- dt_all[dt_all$knockout %in% c(chet, homs)]
    dt_chet <- dt_all[dt_all$knockout == chet,]
    dt_hom <- dt_all[dt_all$knockout == homs,]

    dt_all_per_sample <- aggregate(pKO ~ s, data = dt_all, FUN = length)
    colnames(dt_all_per_sample)[2] <- "all"
    dt_ko_per_sample <- aggregate(pKO ~ s, data = dt_ko, FUN = length)
    colnames(dt_ko_per_sample)[2] <- "ko"
    dt_chet_per_sample <- aggregate(pKO ~ s, data = dt_chet, FUN = length)
    colnames(dt_chet_per_sample)[2] <- "chet"
    dt_hom_per_sample <- aggregate(pKO ~ s, data = dt_hom, FUN = length)
    colnames(dt_hom_per_sample)[2] <- "hom"
    dt_het_per_sample <- aggregate(pKO ~ s, data = dt_het, FUN = length)
    colnames(dt_het_per_sample)[2] <- "het"
    dt_cis_per_sample <- aggregate(pKO ~ s, data = dt_cis, FUN = length)
    colnames(dt_cis_per_sample)[2] <- "cis"

    mymerge = function(x,y) merge.data.table(x,y,all=TRUE)
    mrg <- Reduce(
        mymerge,
           list(
               dt_all_per_sample,
               dt_ko_per_sample, 
               dt_chet_per_sample, 
               dt_hom_per_sample, 
               dt_het_per_sample, 
               dt_cis_per_sample
           )
          )
    mrg[is.na(mrg)] <- 0
    colnames(mrg)[colnames(mrg) == "s"] <- "eid"
    outfile <- paste0(args$out_prefix, ".unrel.counts.txt.gz")
    fwrite(mrg, outfile)

    # merge with phenotypes
    final <- merge(ph, mrg, by = 'eid')
    nrow(final)

    # make sure that encodign is correct
    #str(final[,colnames(final) %in% covars, with = FALSE])
    final$ukbb.centre <- factor(final$ukbb.centre)
    final$sex <- factor(final$sex)

    # run all phenotypes
    phenotypes <- header$V1[header$V1 %in% colnames(ph)]
    phenotypes <- phenotypes[1:300]
    v <- as.character(args$variable)
    print(v)
    print(colnames(final))

    if (!(v %in% colnames(final))) stop(paste(v, "is not among column names!"))
    fits <- lapply(phenotypes, function(pheno){
        write(pheno, stderr())
        # setup model
        if (is.logical(final[[pheno]])){
            covariates <- paste0(covars, collapse = "+")
            model_str <- paste0(pheno ,"~",v,"+", covariates)
            model <- as.formula(model_str)
            # run model
            fit <- glm(
                formula = model,
                data = final,
                family = binomial(link="logit")
            )
            return(fit)
        }
    })

    # setup names
    names(fits) <- phenotypes
    coeffecients <- do.call(rbind, lapply(names(fits), function(fit_name){
        f <- fits[[fit_name]]
        ko_coef <- t(coef(summary(f))[2,])
        ko_coef <- data.table(ko_coef)
        ko_coef$phenotype <- fit_name
        return(ko_coef)
    }))

    outfile = paste0(args$out_prefix, ".unrel.glm.txt.gz")
    write(paste0("writing to ", outfile), stderr())
    fwrite( coeffecients, outfile, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotypes_path", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--phenotypes_header", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--unrelated_path", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--covars_path", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--knockout_dir", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--knockout_pattern", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--variable", default="ko", required = TRUE, help = "what variable to regress (ko, chet, het, homs)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


