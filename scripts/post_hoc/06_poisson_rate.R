
library(data.table)
library(argparse)
library(ggplot2)
library(MASS)

main <- function(args){
  
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
    models <- c("is_chet" , "is_hom", "is_ko", "is_het","is_cis")
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
          if (anno == "damaging_missense") model <-  paste0(m,"~x+mis+frameshift")
          if (anno == "pLoF") model <-  paste0(m,"~x+frameshift+splice_site+non")
          if (anno == "pLoF_damaging_missense") model <-  paste0(m,"~x+mis+frameshift+splice_site+non")
          
          lst[[m]][[gene_set]][[anno]] <- list()
          print(anno)
          print(model)
          # setup fit
          dt_fit <- dt[dt$annotation %in% anno,]
          background_fit <- background[!background$gene_id %in% dt_fit$gene_id,]
          dt_fit <- rbind(dt_fit, background_fit)
          # add background
          dt_fit$x <- FALSE
          dt_fit$x[dt_fit$gene_id %in% genes_in_geneset ] <- TRUE
          # fit models
          fit <- glm(as.formula(model), data = dt_fit, family=poisson(link="log"), control = glm.control(maxit = 200))
          fit_coef <- data.frame(coef(summary(fit)))
          colnames(fit_coef) <- c("est", "error", "z", "p")
          conf <- suppressMessages(exp(cbind(coef(fit), confint(fit))))
          # setup coeffecients
          fit_coef$ci_est <- conf[,1]
          fit_coef$ci_lower <- conf[,2]
          fit_coef$ci_upper <- conf[,3]
          fit_coef$keep <- as.logical(c(0, 1, rep(0, nrow(conf)-2)))
          fit_coef$geneset <- gene_set
          fit_coef$annotation <- anno
          fit_coef$model <- m
          lst[[m]][[gene_set]][[anno]] <- fit_coef
        }
      }
    }


    # combine all of the data.frames
    combined <- do.call(rbind, lapply(lapply(lst, function(l) lapply(l, rbindlist)), rbindlist))
    combined <- combined[combined$keep, ]
    combined$geneset <- factor(combined$geneset, levels=rev(unique(combined$geneset)))
    combined$geneset_cat <- ifelse(combined$geneset %in% essential_names, "Essential", "Non-essential")
    combined <- combined[combined$geneset %in% c(essential_names, non_essential_names),]
    combined$model <- factor(combined$model, levels = models)

    # writ efile
    outfile <- paste0(args$out_prefix,".txt")
    write(paste("writing", outfile), stderr())
    fwrite(combined, outfile, sep = "\t")

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
parser$add_argument("--dir_genesets", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)


