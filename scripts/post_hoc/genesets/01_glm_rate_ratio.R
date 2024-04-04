
library(data.table)
library(argparse)
library(ggplot2)
library(MASS)
library(pscl)

fit_glm_poisson <- function(model, dt_fit, gene_set, anno, m) {
  fit <- glm(as.formula(model), data = dt_fit, family = poisson(link = "log"), control = glm.control(maxit = 200))
  fit_coef <- data.frame(coef(summary(fit)))
  colnames(fit_coef) <- c("est", "error", "z", "p")
  conf <- suppressMessages(exp(cbind(coef(fit), confint(fit))))
  
  fit_coef$ci_lower <- conf[, 1]
  fit_coef$ci_upper <- conf[, 2]
  fit_coef$geneset <- gene_set
  fit_coef$annotation <- anno
  fit_coef$model <- m
  
  return(fit_coef)
}

fit_zinb <- function(model, dt_fit, gene_set, anno, m) {
  zi_model <- zeroinfl(as.formula(paste(model, "| 1")), data = dt_fit, dist = "negbin")
  model_summary <- summary(zi_model)
  
  coef_summary_count <- model_summary$coefficients$count
  exp_conf_ints <- exp(confint(zi_model))
  
  exp_est = exp(coef_summary_count[, "Estimate"])
  fit_coef <- data.frame(
    est = coef_summary_count[, "Estimate"],
    error = coef_summary_count[, "Std. Error"],
    z = coef_summary_count[, "z value"],
    p = coef_summary_count[, "Pr(>|z|)"],
    exp_est = exp_est,
    ci_lower = exp_conf_ints[, 1],
    ci_upper = exp_conf_ints[, 2]
  )
  
  fit_coef$geneset <- gene_set
  fit_coef$annotation <- anno
  fit_coef$model <- m
  
  return(fit_coef)
}


main <- function(args){

    print(args)

    # remove the three outliers
    #aggr_mrg <- fread("~/Downloads/combined_annotations_by_sample.new.counts.txt.gz")
    cols_aggr <- c("gene_id", "transcript_id", "annotation")
    aggr_mrg <- fread(args$in_count)
    genes_to_exclude <- c("ENSG00000094963","ENSG00000178917","ENSG00000188163")
    aggr_mrg <- aggr_mrg[!aggr_mrg$gene_id %in% genes_to_exclude,]

    # create synthetic category
    aggr_sx_chet <- aggregate(is_chet ~ gene_id + transcript_id, data = aggr_mrg, FUN = sum)
    aggr_sx_hom <- aggregate(is_hom ~ gene_id + transcript_id, data = aggr_mrg, FUN = sum)
    aggr_sx_cis <- aggregate(is_cis ~ gene_id + transcript_id, data = aggr_mrg, FUN = sum)
    aggr_sx_het <- aggregate(is_het ~ gene_id + transcript_id, data = aggr_mrg, FUN = sum)
    aggr_sx_mrg <- merge(aggr_sx_chet, aggr_sx_hom, by = cols_aggr[1:2], all = TRUE)
    aggr_sx_mrg <- merge(aggr_sx_mrg, aggr_sx_cis, by = cols_aggr[1:2], all = TRUE)
    aggr_sx_mrg <- merge(aggr_sx_mrg, aggr_sx_het, by = cols_aggr[1:2], all = TRUE)
    aggr_sx_mrg$annotation <- "combined"
    aggr_sx_mrg <- aggr_sx_mrg[,colnames(aggr_mrg)]
    aggr_sx_mrg$is_ko <- aggr_sx_mrg$is_chet + aggr_sx_mrg$is_hom

    # create mapping for knockout burden
    mapping_chet <- aggr_sx_mrg$is_chet
    names(mapping_chet) <- aggr_sx_mrg$gene_id
    mapping_hom <- aggr_sx_mrg$is_hom
    names(mapping_hom) <- aggr_sx_mrg$gene_id
    mapping_ko <- aggr_sx_mrg$is_ko
    names(mapping_ko) <- aggr_sx_mrg$gene_id
    mapping_het <- aggr_sx_mrg$is_het
    names(mapping_het) <- aggr_sx_mrg$gene_id

    # combine
    aggr_mrg <- aggr_mrg
    aggr_mrg$cond_chet <- mapping_chet[aggr_mrg$gene_id]
    aggr_mrg$cond_hom <- mapping_hom[aggr_mrg$gene_id]
    aggr_mrg$cond_ko <- mapping_ko[aggr_mrg$gene_id]
    aggr_mrg$cond_het <- mapping_het[aggr_mrg$gene_id]
    aggr_mrg$total <- aggr_mrg$is_chet + aggr_mrg$is_hom + aggr_mrg$is_cis
    aggr_mrg$is_ko <- aggr_mrg$is_chet + aggr_mrg$is_hom

    # get list of background genes
    background <- aggr_mrg
    background$annotation <- NA
    background$is_chet <- 0
    background$is_hom <- 0
    background$is_cis <- 0
    background$is_ko <- 0
    background$is_het <- 0
    background$total <- 0
    background <- background[!duplicated(background),]

    # get mutation rates (note some genes are not defined)
    #mr <- fread("~/Projects/09_genesets/genesets/data/mutation_rates/samocha2014.txt.gz")
    mr <- fread(args$file_mutation_rates)
    colnames(mr)[colnames(mr) == "ensembl_gene_id"] <- "gene_id"
    mr$transcript <- NULL
    mr$hgnc_symbol <- NULL
    dt <- merge(aggr_mrg, mr, all.x =TRUE)
    background <- merge(background, mr, all.x = TRUE)

    # get essential non-essential genes
    files <- list.files(args$dir_genesets, full.names = TRUE)
    files <- files[!grepl(files, pattern="essential_non_essential_combined")]
    gene_lst <- lapply(files, function(x){fread(x)$gene_id})
    names(gene_lst) <- tools::file_path_sans_ext(basename(files))
        
    essential_names <- c(
      "essential_in_mice_georgi2013",
      "essential_gnomad_karczewski2020",
      "essential_adam_vincenti2021",
      "essential_in_culture_hart2014",
      "essential_crispr_hart2017",
      "pLI>=0.9"
      )

    non_essential_names <- c(
      "non_essential_gnomad_karczewski2020",
      "non_essential_in_culture_hart2014",
      "homozygous_lof_tolerant_karczewski2020"
    )

    # get loeuf
    #loeuf <- fread("~/Projects/09_genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv")
    loeuf <- fread(args$file_pli)
    gene_lst[["pLI>=0.9"]] <- unique(as.character(na.omit(loeuf$gene_id[loeuf$pLI >= 0.9])))

    # omim gene list
    if (!is.null(args$file_omim)){
        #omim <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/omim/230329_morbidmap_by_gene_with_inheritance.txt")
        omim <- fread(args$file_omim)
        omim <- omim[grepl(omim$ensgid, pattern="ENSG"),]
        omim_autosomal_recessive <- unique(omim$ensgid[omim$AR])
        omim_autosomal_dominant <- unique(omim$ensgid[omim$AD])
        gene_lst[["omim_AR"]] <- omim_autosomal_recessive
        gene_lst[["omim_AD"]] <- omim_autosomal_dominant
    } 
    
    # GTEx files
    if (!is.null(args$file_gtex)){
        
        # load gtex and get top 10% specific
        #d <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gtex/GTEx.tstat.tsv")
        d <- fread(args$file_gtex)
        d <- cbind(gene_id = d$ENSGID, data.table(apply(d[,-1], 2, function(x) x > quantile(x, probs = 0.90))))
        colnames(d) <- paste0("gtex_",tolower(gsub("\\_$","",gsub("\\_+","\\_",gsub("(\\()|(\\)|(\\-))","_",colnames(d))))))

        # add to gene_lst 
        tissues <- colnames(d)[-1]
        for (tissue in tissues){
            bool_tissue <- d[[tissue]]
            genes_specific_to_tissue <- d$gtex_gene_id[bool_tissue]
            gene_lst[[tissue]] <- genes_specific_to_tissue
        }
    }

    # cancer gene set
    if (!is.null(args$file_cancer)){
        d <- fread(args$file_cancer)
        #d <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gsea/cancer/gavish_3ca_genes.txt")
        d <- d[(!is.na(d$ensgid)) & (d$ensgid != "")]
        genesets <- unique(d$pathway)
        for (geneset in genesets[1]){
            genes <- d$ensgid[d$pathway == geneset]
            gene_lst[[geneset]] <- genes
        }
    }

    #models <- c("is_chet" , "is_hom", "is_ko", "is_het")
    #models <- c("is_cis")
    valid_models <- c("is_chet" , "is_hom", "is_ko", "is_het","is_cis")
    models <- unlist(strsplit(args$models_to_run, split=","))
    stopifnot(all(models %in% valid_models))
    
    annotations <- unique(aggr_mrg$annotation)

    get_covars <- function(model){
      return(unlist(strsplit(model, split = "\\+"))[-1])
    }

    lst <- list()
    for (m in models){
      lst[[m]] <- list()
      model <- paste0(m,"~x+bp")
      for (gene_set in names(gene_lst)){
        lst[[m]][[gene_set]] <- list()
        genes_in_geneset <- gene_lst[[gene_set]]
        print(gene_set)
        for (anno in annotations){
          
          if (anno == "synonymous") model <- paste0(m,"~x+syn")
          if (anno == "other_missense") model <-  paste0(m,"~x+mis")
          if (anno == "damaging_missense") model <-  paste0(m,"~x+mis")
          if (anno == "pLoF") model <-  paste0(m,"~x+frameshift+splice_site+non")
          if (anno == "pLoF_damaging_missense") model <-  paste0(m,"~x+mis+frameshift+splice_site+non")
          
          #print("USING ALL MODEL!!!")
          #model <- paste(m,"~x+all")

          lst[[m]][[gene_set]][[anno]] <- list()
          print(anno); print(model)
          
          # setup fit
          dt_fit <- dt[dt$annotation %in% anno,]
          background_fit <- background[!background$gene_id %in% dt_fit$gene_id,]
          dt_fit <- rbind(dt_fit, background_fit)
          
          # add background
          dt_fit$x <- FALSE
          dt_fit$x[dt_fit$gene_id %in% genes_in_geneset ] <- TRUE

          # trasnforming covariates seem to solve convergence issue
          # see https://stats.stackexchange.com/questions/76488/error-system-is-computationally-singular-when-running-a-glm         
          scaling <- 100
          dt_fit$syn <- dt_fit$syn * scaling
          dt_fit$mis <- dt_fit$mis * scaling
          dt_fit$frameshift <- dt_fit$frameshift * scaling
          dt_fit$splice_site <- dt_fit$splice_site * scaling
          dt_fit$non <- dt_fit$non * scaling

          # fit GLM
          if (args$glm_method == "poisson") {
                fit_coef <- fit_glm(model, dt_fit, gene_set, anno, m)
          } else if (args$glm_method == "zinfb") {
                fit_coef <- fit_zinb(model, dt_fit, gene_set, anno, m)
          } else {
                stop(paste(args$glm_method), "is not a valid method to use")
          }

          # only keep variable of interest
          fit_coef$keep <- rownames(fit_coef) %in% "xTRUE"
          lst[[m]][[gene_set]][[anno]] <- fit_coef

        }
      }
    }

    # combine all of the data.frames
    combined <- do.call(rbind, lapply(lapply(lst, function(l) lapply(l, rbindlist)), rbindlist))
    combined <- combined[combined$keep, ]
    combined$geneset <- factor(combined$geneset, levels=rev(unique(combined$geneset)))
    #combined$geneset_cat <- ifelse(combined$geneset %in% essential_names, "Essential", "Non-essential")
    #combined <- combined[combined$geneset %in% c(essential_names, non_essential_names),]
    combined$model <- factor(combined$model, levels = models)

    # write efile
    outfile <- paste0(args$out_prefix,".txt")
    write(paste("writing", outfile), stderr())
    fwrite(combined, outfile, sep = "\t", na="NA")

    # setup colors
    categories <- c('pLoF','pLoF_damaging_missense','damaging_missense','other_missense','synonymous', 'combined')
    my_colors <- c("#B13F64","#DD686D","#F09D7C", "#F4D400", "#7CA98A", "grey30")
    names(my_colors) <- categories
    fill_scale <- scale_fill_manual(name = "annotation", values = my_colors)
    color_scale <- scale_color_manual(name = "annotation", values = my_colors)
    combined$annotation <- factor(combined$annotation, levels = categories)


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--in_count", default=NULL, required = TRUE, help = "")
parser$add_argument("--file_mutation_rates", default=NULL, required = TRUE, help = "")
parser$add_argument("--file_pli", default=NULL, required = TRUE, help = "")
parser$add_argument("--file_omim", default=NULL, required = FALSE, help = "")
parser$add_argument("--file_cancer", default=NULL, required = FALSE, help = "")
parser$add_argument("--file_gtex", default=NULL, required = FALSE, help = "")
parser$add_argument("--dir_genesets", default=NULL, required = TRUE, help = "")
parser$add_argument("--glm_method", default="zinfb", required = TRUE, help = "either 'zinfb' or 'poisson'")
parser$add_argument("--models_to_run", default="is_chet", required = TRUE, help = "either 'is_hom', 'is_chet', 'is_het', 'is_cis' or 'is_ko'")
args <- parser$parse_args()

main(args)


