#!/usr/bin/env Rscript

devtools::load_all('utils/modules/R/phasingtools')
library(argparse)
library(RColorBrewer)
library(stringr)

main <- function(args){

    # parser
    print(args)
    stopifnot(dir.exists(args$in_dir_merged))     

    # load files in directories
    dirs <- strsplit(args$in_dir, split = ',')
    stops <- lapply(dirs, function(d) stopifnot(dir.exists(d)))
    files <- unlist(lapply(dirs, function(d){list.files(d, pattern = 'trio', full.names = TRUE)}))
    
    # organize by filename, chunks and extension 
    M <- as.data.frame(do.call(rbind, strsplit(basename(files), split = '\\.')))
    colnames(M) <- c('param','n','ext')
    datasets = unique(M$param)

    # read each file and calculate average contribution of chunk
    comparison <- do.call(rbind, lapply(datasets, function(dname){
    
      regex_bool <- paste0('^',dname)
      selected_files <- files[grepl(regex_bool, basename(files))]
       
      ds <- lapply(selected_files, function(f) {
            d <- suppressWarnings(summarize_bcftools_trio_stats(fread(f)))
            d$chunk <- unlist(strsplit(str_match(f, ".of."), split = 'of'))[1]
            d$dataset <- dname
            return(d)
        })
        
      ds <- as.data.frame(do.call(rbind, ds))
      
      avg <- tail(ds)[1,]
      avg$n_tested <- sum(ds$n_tested)
      avg$n_switch <- sum(ds$n_switch)
      avg$n_mendel <- sum(ds$n_mendel)
      bconf <- Hmisc::binconf(avg$n_switch, avg$n_tested)
      est <- bconf[1]
      err <- abs(est - bconf[3])
      avg$ci_ser_est <- est
      avg$ci_ser_error <- err
      avg$ci_ser_est_pct <- NA
      avg$ci_ser_error_pct <- NA
      avg$chunk <- as.character('Average')
      
      ds <- rbind(ds, avg)
      
      return(ds)
    
    }))

    # get results from merged files
    # works only for chr20 for now! 
    files_merge = list.files(args$in_dir_merged, 
                   pattern = 'trio', 
                   full.names = TRUE) 
    
    dfs_merge <- do.call(rbind, lapply(files_merge, function(f){
        stopifnot(file.exists(f))
        merged <- suppressWarnings(summarize_bcftools_trio_stats(fread(f)))
        merged$dataset <- gsub("_20", "", tools::file_path_sans_ext(basename(f)))
        merged$chunk <- as.character("Merged")
        return(merged)
    }))

    # combine data
    combined <- rbind(comparison, dfs_merge)

    # assume maximum of 4 chunks (index 5 is avg, and 6 is merge)
    myColors <- brewer.pal(6,"Set1")
    myColors[5] <- 'grey'
    myColors[6] <- 'black'
    names(myColors) <- levels(combined$chunk)
    colScale <- scale_colour_manual(name = "chunk",values = myColors)

    # get values for vertical lines
    best_mrg <- combined[combined$chunk == 'Merged',]
    best_mrg <- best_mrg[best_mrg$ci_ser_est == min(best_mrg$ci_ser_est),]
    best_avg <- combined[combined$chunk == 'Average',]
    best_avg <- best_avg[best_avg$ci_ser_est == min(best_avg$ci_ser_est),]
   
    # convert into ordered factor (for plotting)
    combined$chunk <- factor(combined$chunk, levels = c("1","2","3","4","Average","Merged"))
    combined[is.na(combined$chunk),]

    # plot resulst
    pd <- position_dodge(0.7)
    plt <- ggplot(combined, 
           aes(
               x=100*ci_ser_est, 
               xmax = 100*(ci_ser_est + ci_ser_error),
               xmin = 100*(ci_ser_est - ci_ser_error),
               y = dataset,
               color = chunk 
           )) +
        colScale + 
        geom_vline(xintercept=best_mrg$ci_ser_est*100, linetype = 'dashed', col = 'black') +
        geom_vline(xintercept=best_avg$ci_ser_est*100, linetype = 'dashed', color = 'black') +
        geom_point(stat='identity', position = pd, size = 2) +
        geom_errorbar(stat='identity', position = pd,width = 0.75) +
        ggtitle('Tradeoff between phasing region size and phasing accuracy',) +
        labs(color = "Phasing chunk(s)") +
        xlab('Switch Errors (%)') + ylab('') + 
        theme_bw()

    # write out plot and table
    outfile_plt = paste0(args$out_prefix, ".pdf")
    outfile_tbl = paste0(args$out_prefix, ".tsv")
    fwrite(combined, outfile_tbl, sep = '\t')
    pdf(outfile_plt, width = 8, height = 6)
    print(plt)
    graphics.off()
    #ggplot(filename = outfile_plt, plot = plt, width = 8, height = 6)
   
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--in_dir_merged", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

