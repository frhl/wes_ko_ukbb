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
    dirs <- unlist(strsplit(args$in_dir, split = ','))
    dirs <- dirs[dirs != ""]
    #dirs[[1]] <- dirs[[1]][dirs[[1]] != ""]
    #print(dirs)
    stops <- lapply(dirs, function(dir) {
        if (!dir.exists(dir)) stop(paste0(dir, ' (directory) does not exist!'))
    })
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
            d$chr <- str_extract(f, 'chr[0-9]+')
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
        
        merged$dataset <- tools::file_path_sans_ext(basename(f))
        merged$chunk <- as.character("Merged")
        merged$chr <- str_extract(f, 'chr[0-9]+')
        return(merged)
    }))

    combined <- rbind(comparison, dfs_merge)

    max_chunks <- max(na.omit(as.numeric(as.character(combined$chunk))))
    myColors <- brewer.pal(max_chunks,"PiYG")
    myColors[10] <- 'grey'
    myColors[11] <- 'black'
    names(myColors) <- levels(combined$chunk)
    colScale <- scale_colour_manual(name = "chunk",values = myColors)

    # get values for vertical lines
    best_mrg <- combined[combined$chunk == 'Merged',]
    best_mrg <- best_mrg[best_mrg$ci_ser_est == min(best_mrg$ci_ser_est),]
    #best_avg <- combined[combined$chunk == 'Average',]
    #best_avg <- best_avg[best_avg$ci_ser_est == min(best_avg$ci_ser_est),]

    # convert into ordered factor (for plotting)
    combined$chunk <- factor(combined$chunk, levels = c(as.character(1:9),"Average","Merged"))
    combined$chr <- factor(combined$chr, levels = paste0("chr",1:22))

    # plot resulst
    pd <- position_dodge2(width = 0.6, preserve = "single")
    ggplot(combined,
           aes(
               y=100*ci_ser_est,
               ymax = 100*(ci_ser_est + ci_ser_error),
               ymin = 100*(ci_ser_est - ci_ser_error),
               x = chunk,
               fill = chunk
           )) +
        fillScale +
        geom_bar(stat='identity', color = 'black') +
        geom_hline(yintercept=best_lig$ci_ser_est*100, linetype = 'dashed', col = 'black', alpha = 0.8) +
        geom_hline(yintercept=worst_lig$ci_ser_est*100, linetype = 'dashed', col = 'black', alpha = 0.8) +
        geom_errorbar(stat='identity', position = pd, width = 0.4) +
        ggtitle('Phasing accuracy across chunks (UKBB Calls+WES)',
                'Dashed line indicate lowest and highest SERs observed. M = Merged, L = Ligated.') +
        labs(color = "Phasing chunk(s)") +
        ylab('Switch Errors (%)') + xlab('Chromosmal chunks') +
        theme_bw() +
        facet_wrap(~chr) 
    
    #pd <- position_dodge(0.7)
    #plt <- ggplot(combined,
    #       aes(
    #           x=100*ci_ser_est,
    #           xmax = 100*(ci_ser_est + ci_ser_error),
    #           xmin = 100*(ci_ser_est - ci_ser_error),
    #           y = chr,
    #           color = chunk
    #       )) +
    #    colScale +
    #    geom_vline(xintercept=best_mrg$ci_ser_est*100, linetype = 'dashed', col = 'black') +
    #    #geom_vline(xintercept=best_avg$ci_ser_est*100, linetype = 'dashed', color = 'black') +
    #    geom_point(stat='identity', position = pd, size = 2) +
    #    geom_errorbar(stat='identity', position = pd,width = 0.75) +
    #    ggtitle('Phasing accuracy across UKBB WES+CALLS',) +
    #    labs(color = "Phasing chunk(s)") +
    #    xlab('Switch Errors (%)') + ylab('') +
    #    theme_bw()

    # write out plot and table
    outfile_plt = paste0(args$out_prefix, ".pdf")
    outfile_tbl = paste0(args$out_prefix, ".tsv")
    fwrite(combined, outfile_tbl, sep = '\t')
    pdf(outfile_plt, width = args$img_width, height = args$img_height)
    print(plt)
    graphics.off()
    
    # write error rates by site
    outfile_sites <- paste0(args$out_prefix, "sites.pdf")
    pdf(outfile_sites, width = 10, height = 8)
    for (chr in 1:22){
        #chunk_dir <- "data/phased/wes_union_calls/chunks/final/ukb_eur_wes_union_calls_200k_chrCHR-16xshort.qe"
        chunk_dir <- gsub("CHR",chr, args$chunk_dir)
        chunks <- read_chunks_from_dir(chunk_dir)
        
        ligated_files <- list.files(args$ligated_dir, full.names = TRUE, pattern = '[0-9]+.txt')
        ligated_file <- ligated_files[grepl(paste0("chr",chr,"-"), ligated_files)]
        ligated <- read_chunks_combined(ligated_file)
        ligated$chunk <- "ligated"
        
        merged_files <- list.files(args$merged_dir, full.names = TRUE, pattern = '[0-9]+.txt')
        merged_file <- merged_files[grepl(paste0("chr",chr,"-"), merged_files)]
        merged <- read_chunks_combined(merged_file)
        merged$chunk <- "ligated"

        mrg <- rbind(chunks, ligated, merged)
        max_chunks <- suppressWarnings(max(na.omit(as.numeric(as.character(mrg$chunk)))))
        mrg$chunk <- factor(mrg$chunk, levels = as.character(c(1:max_chunks, "merged","ligated")))
   
        #trim <- fread('data/phased/wes_union_calls/trimmed/ukb_eur_wes_union_calls_200k_chr10_trims.txt')
        trim_files <- list.files(args$trim_dir, full.names = TRUE, pattern = '[0-9]+.txt')
        trim_file <- trim_files[grepl(paste0("chr",chr,"_"), trim_files)]
        trim <- fread(trim_file)
        colnames(trim) <- c('dir','prefix','i','mt1_right_flank',  'mt2_left_flank','o','start','end')
        trim$midpoint <- trim$mt1_right_flank-trim$mt2_left_flank 

        colors <- brewer.pal(max_chunks,"Blues")
        colors[max_chunks + 1] <- 'orange'
        colors[max_chunks + 2] <- 'darkred'
        names(colors) <- levels(mrg$chunk)
        scale_color <- scale_colour_manual(name = "chunk",values = colors)

        plt <- ggplot(mrg, aes(x = POS, y = cumsum, group = chunk, color = chunk)) +
            geom_point(alpha = 0.8) +
            geom_line() +
            scale_color +
            annotate("rect", xmin = trim$mt2_left_flank, xmax = trim$mt1_right_flank, 
                     ymin = -Inf, ymax = Inf, alpha = .2) +
            geom_vline(xintercept=trim$mt1_right_flank, linetype = 'dashed', color = 'black') +
            geom_vline(xintercept=trim$mt2_left_flank, linetype = 'dashed', color = 'black') +
            xlab("BP (genomic position)") +
            ylab("Cumulativ sum of Switch Errors") +
            theme_bw() +
            ggtitle(paste("Switch Errrors when ligating Chromosome", chr),
                    "Shaded regions indicate ligation variant overlap")

        print(plt)
     }
    graphics.off()   


#ggplot(filename = outfile_plt, plot = plt, width = 8, height = 6)
   
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--ligated_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--merged_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--chunk_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--trim_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--in_dir_merged", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--img_width", default=8, help = "Where should the results be written?")
parser$add_argument("--img_height", default=6, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

