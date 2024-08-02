
library(data.table)
library(argparse)
library(ggplot2)
library(scales)

# get table of variants by pp to plot
get_seq_table <- function(data){
    xs <- seq(0.50, 0.999, by = 0.001)
    ys <- unlist(lapply(xs, function(x) return(sum(data$PP > x))))
    table <- data.table(pp_cutoff = xs, retained = ys)
    solver <- function(pp) return(sum(data$PP > pp))
    return(list(table = table, solve = solver))
}


# actually plot variants by pp
ggplot_ppxn <- function(dt, main ){
    # get table of things to plot
    lst <- get_seq_table(dt)
    table <- lst$table
    f <- lst$solve
    # actually plot it
    p <- ggplot(table, aes(x=pp_cutoff, y=retained)) +
        geom_line() +
        geom_vline(xintercept = 0.99, col = "blue", linetype = 'dashed') +
        geom_hline(yintercept = f(0.99), col = "blue", linetype = 'dashed') +
        geom_vline(xintercept = 0.90, col = "red", linetype = 'dashed') +
        geom_hline(yintercept = f(0.90), col = "red", linetype = 'dashed') +
        ylab("log10(Variants retained)") +
        xlab("PP Cutoff") +
        scale_y_continuous(trans = 'log10',
                           breaks = trans_breaks('log10', function(x) 10^x),
                           labels = trans_format('log10', math_format(10^.x))) +
        theme_bw() +
        ggtitle(main) +
        theme(
            strip.text = element_text(size=14),
            axis.text=element_text(size=12),
            axis.title=element_text(size=12,face="bold"),
            plot.title = element_text(hjust=0.5),
            plot.subtitle = element_text(hjust=0.5),
            legend.position="none"
        ) 
    return(list(plt = p, data = table))
}



main <- function(args){
   
    paths <- args$in_dir
    pattern <- args$regex
    out_prefix <- args$out_prefix

    files <- list.files(paths, pattern = pattern, full.names = TRUE)
    print(files)
    stopifnot(length(files)>21)
    dt <- do.call(rbind,lapply(files, function(f){
        msg <- paste0("Reading ", f)
        write(msg, stderr())
        d <- fread(f)
        d$AC <- as.numeric(gsub("(\\[)|(\\])", "", d$AC))
        d$max_AN <- max(d$AC)
        d$singleton <- (d$AC == 1) | (d$AC == d$max_AN)
        d$doubleton <- (d$AC == 2) | (d$AC == (d$max_AN-1))
        d$tripleton <- (d$AC == 3) | (d$AC == (d$max_AN-2))
        return(d)
    }))

    # plot and get data that is beign plotted
    lst1 <- ggplot_ppxn(dt[dt$singleton == TRUE], "Singletons")
    lst2 <- ggplot_ppxn(dt[dt$doubleton == TRUE], "Doubletons")
    lst3 <- ggplot_ppxn(dt[dt$tripleton == TRUE], "Tripletons")    
    lst4 <- ggplot_ppxn(dt, "All")    

    outfile1 <- paste0(out_prefix, "_singletons")
    ggsave(paste(outfile1,".png"), lst1$plt, width = 8, height = 8)
    fwrite(lst1$data, paste(outfile1, ".txt.gz"))

    outfile2 <- paste0(out_prefix, "_doubletons")
    ggsave(paste(outfile2,".png"), lst2$plt, width = 8, height = 8)
    fwrite(lst2$data, paste(outfile2, ".txt.gz"))

    outfile3 <- paste0(out_prefix, "_tripletons")
    ggsave(paste(outfile3,".png"), lst3$plt, width = 8, height = 8)
    fwrite(lst3$data, paste(outfile3, ".txt.gz"))

    outfile4 <- paste0(out_prefix, "_all")
    ggsave(paste(outfile4,".png"), lst4$plt, width = 8, height = 8)
    fwrite(lst4$data, paste(outfile4, ".txt.gz"))





}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--regex", default=NULL, required = TRUE, help = "Chromosome")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


