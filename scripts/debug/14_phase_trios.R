library(data.table)
library(Hmisc)
setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')
devtools::load_all('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/ukbtools')

# setup parser
parser <- ArgumentParser()
parser$add_argument("--in_file", default=NULL, help = "the path to the input file")
parser$add_argument("--out_file", default=NULL, help = "the path to the output")
args <- parser$parse_args()

# load in main file and calculate binomial CI
d <- zcat(args$in_file)
n <- nrows(d)
alpha = 0.05
lst < lapply(1:n, function(i){
       row <- d[i, ,with = FALSE]
       bconf <- as.data.table(Hmisc::bincomf(row$errors_sum, row$hetz_sum, alpha = alpha))
       bconf$alpha <- alpha
       return(bconf)
})

# writingn to out file
d_out <- cbind(d, do.call(rbind, lst))
fwrite(d_out, args$out_file, row.names = FALSE, quote = FALSE)
write(paste('writing to', args$out_file),stdout())



