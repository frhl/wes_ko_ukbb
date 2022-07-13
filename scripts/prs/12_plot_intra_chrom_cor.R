
library(argparse)
library(data.table)
library(ggplot2)

discard_ext <- function(x) tools::file_path_sans_ext(x)
label_sig <- function(x, p){return(ifelse(x < 0.05, ifelse(x < p, "**","*"),""))}

main <- function(args){

    #files <- list.files("data/prs/validation/chrom/old", full.names = TRUE)
    files <- list.files(args$in_Dir, full.names = TRUE)
    lst <- lapply(files, function(f){
        d <- fread(f)
        d$path <- gsub("_pgs_chrom","",unlist(strsplit(basename(f), '\\.'))[1])
        dup <- d[,c("col1","col2")]
        dups <- duplicated(do.call(rbind, lapply(1:nrow(dup), function(i) sort(as.numeric(gsub("chr","",dup[i]))))))
        d <- d[!dups,]
        return(d)
    })  

    # combine list of data.table
    d <- do.call(rbind, lst)
    autosomes <- paste0("chr",1:22)
    d$label <- label_sig(d$pvalue, 0.05/231)
    # plot them
    plt <- ggplot(d, 
           aes(
               x=factor(col1, levels = autosomes), 
               y=factor(col2, levels = autosomes), 
               label = label,
               fill=correlation +
        geom_tile() + 
        geom_text() + 
        theme_bw() + 
        xlab("") +
        ylab("") +
        #labs(fill=expression(-log[10]~"(P-value)")) +
        labs(fill="Pearson correlation") +
        scale_fill_gradient2(mid="white",low = "blue", high = 'red', limits = c(-0.1, 0.1)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        theme(aspect.ratio = 1) +
        facet_wrap(~path, scales = "free")

    outfile = paste0(args$out_prefix,".png")
    ggsave(plt, outfile, width = args$width, height = args$height)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "in directory")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--out_height", default=NULL, required = TRUE, help = "image height")
parser$add_argument("--out_width", default=NULL, required = TRUE, help = "image width")
args <- parser$parse_args()

main(args)









