
library(data.table)
library(ggplot2)


path <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/knockouts/tables/poisson_rate_essential_genesets_review.txt"
combined <- fread(path)

# setup colors
categories <- c('pLoF','pLoF_damaging_missense','damaging_missense','other_missense','synonymous', 'non_coding')
my_colors <- c("#B13F64","#DD686D","#F09D7C", "#F4D400", "#7CA98A", "#C2D4D8")
combined <- combined[combined$annotation %in% categories,]
names(my_colors) <- categories
fill_scale <- scale_fill_manual(name = "annotation", values = my_colors)
color_scale <- scale_color_manual(name = "annotation", values = my_colors)
combined$annotation <- factor(combined$annotation, levels = categories)

## mapping for pretty plots

# map back and forth to categories
simplify_gtex_name <- function(x) {gsub("\\.","",gsub(" ","",paste0("gtex_",tolower(gsub("\\_$","",gsub("\\_+","\\_",gsub("(\\()|(\\)|(\\-))","_",x)))))))}
gtex_mapping <- colnames(fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gtex/GTEx.tstat.tsv"))
names(gtex_mapping) <- simplify_gtex_name(gtex_mapping)

# gtex categories 
gtex_cat <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gtex/GTEX.tstat.categories.genoppi.csv")
gtex_cat$Tissue.genoppi <- simplify_gtex_name(gtex_cat$Tissue.genoppi)
# a bit of manual mapping too..
gtex_cat$Tissue.genoppi[gtex_cat$Tissue.genoppi == "gtex_cells_ebvtransformed_lymphocytes"] <- "gtex_cells_ebv_transformed_lymphocytes"
gtex_cat$Tissue.genoppi[gtex_cat$Tissue.genoppi == "gtex_brain_spinal_cord_cervical_c1"] <- "gtex_brain_spinal_cord_cervical_c_1"
gtex_cat_map <- gtex_cat$Tissue.category.for.display
names(gtex_cat_map) <- gtex_cat$Tissue.genoppi

# combine them...
combined_gtex <- combined[grepl(combined$geneset, pattern="gtex"),]
combined_gtex$geneset_original_title <- gtex_mapping[combined_gtex$geneset]
stopifnot(sum(is.na(combined_gtex$geneset_original_title)) == 0)
combined_gtex$geneset_original_category <- gtex_cat_map[combined_gtex$geneset]
stopifnot(sum(is.na(combined_gtex$geneset_original_category)) == 0)

# make the labels prettier
combined_gtex$geneset_category <- combined_gtex$geneset_original_category

# ensure that genestes are orderd
combined_gtex$geneset_label <- gsub("gtex_","", combined_gtex$geneset)
geneset_levels <- unique(combined_gtex$geneset_label[order(combined_gtex$geneset_original_category)])
combined_gtex$geneset_label <- factor(combined_gtex$geneset_label, levels = geneset_levels)

# go over each model and save it
out_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/knockouts/tables"
out_prefix=paste0(out_dir, "/poisson_rate_essential_genesetes_review.gtex")
models <- unique(combined_gtex$model)
pd <- position_dodge(0.7)
for (cur_model in models){
    
    #ggplot(combined_gtex[combined_gtex$annotation %in% "pLoF_damaging_missense",],
    ggplot(combined_gtex[combined_gtex$model %in% "is_hom",],
            aes(
             x=geneset_label,
             y=ci_est,
             ymax=ci_upper,
             ymin=ci_lower,
             fill=geneset_original_category,
             group=annotation
           )
    ) +
      geom_bar(stat="identity", width = 0.5) +
      geom_pointrange(position = pd, size = 0.6) +
      geom_hline(yintercept = 1, linetype = 'dashed', color = "black") +
      scale_y_continuous(
        trans="log10",
        breaks=c(0.01, 0.1, 0.25, 0.5, 1, 3)
      ) +
      #xlim(c(0.001,10)) +
      ylab("Rate ratio") +
      xlab("Geneset") +
      #color_scale +
      theme_bw() +
      labs(shape="Knockout", color = "Annotation") +
      ggtitle(paste0("Poisson regression for ",cur_model," variants")) +
      theme(
        legend.position = "right",
        strip.text = element_text(size=10),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        axis.title.y = element_text(margin=ggplot2::margin(t=10)),
        axis.title.x = element_text(margin=ggplot2::margin(r=10)),
        axis.text.x = element_text(angle=75, hjust=1, size=10),
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5)
      )  +
      facet_grid(annotation~1)
    out_file <- paste0(out_prefix,".", cur_model,".pdf") 
    write(paste("Saving to", out_file), stdout())
    ggsave(out_file,  width = 18, height=14)
}







