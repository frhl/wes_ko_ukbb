# module purge
# conda activate rpy

# setup paths and libs
library(argparse)
library(data.table)
library(ggplot2)
library(dplyr)

#setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, help = "What directory should files be loaded from?")
parser$add_argument("--in_pattern", default = 'genebass', help = "what string should the file contains?")
parser$add_argument("--out_prefix", default = 'plot', help = "string prefix for resulting plots")
args <- parser$parse_args()


# read in .bgz files 
zcat_fread <- function(path,...){
	cmd = paste("zcat", path)
	return(fread(cmd = cmd,...))
}


# expand hail list items. E.g. an array list A = [x,y]
# will be converted into seperate columns, A1 = x and A2 = y.
# returns the full data.frame with the new collumns appended
expand_hail_list <- function(dt, name){
    stopifnot(name %in% colnames(dt))
    x <- dt[[name]]
    x <- gsub('(\\[)|(\\])|','',x)
    mat <- do.call(rbind, strsplit(x, split = ','))
    n <- ncol(mat)
    names <- paste0(name,1:n)
    colnames(mat) <- names
    return(cbind(dt, mat))
}


files = sort(list.files(args$in_dir, pattern = args$in_pattern, full.names = TRUE))
stopifnot(length(files) > 0)
print(files)

pdf(paste0(args$out_prefix,'.pdf'), width = 9, height = 7)
for (f in files){
    print(paste('running',f))
    d = zcat_fread(f)
    d = expand_hail_list(d, name = 'variant_qc.AF')
    p <- ggplot(dt, aes(x=variant_qc.AF2, y=gb_AF)) + 
    geom_bin2d(aes(fill=log10(..count..)), bins=40) + theme_classic() + 
    xlab("UK Biobank AF") + ylab("genebass AF") +
    ggtitle("UK Biobank pan-UKB AF vs. gnomAD AF after MAF matching", basename(f))
    print(p)
    p <- p + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')
    print(p)
}
graphics.off()



