# this script identifies common pLoFs, which we can exclude. The 
# assumption is that if they are this common, they are probably not damaging.

library(fread)

# process variants on each strand
clean_hail_list <- function(x, comma_sub = ', ') {
    x <- gsub('(\\[)|(\\])|(\\")','',x)
    x <- gsub('\\,',comma_sub,x)
    return(x)
}


main <- function(args){
    
    # select only files that pass a threshold
    # note, that we substitute "CHR" when reading in the file
    d <- do.call(rbind, lapply(1:22, function(chr){
        p <- gsub("CHR",chr, args$in_file)
        dt <- fread(p)
        if (nrow(dt) > args$threshold) return(dt)
    }))

    # get the most abundant gene-variant pairs
    table <- data.frame(table(d$gene_id, d$varid))
    table <- table[rev(order(table$Freq)),]
    table <- table[table$Freq > args$threshold,]
    print(paste(nrow(table), 'gene-variants have more than 100k occurences!'))
    colnames(table) <- c('gene_id', 'varid', 'frequency')
    table$varid <- clean_hail_list(table$varid)

    # write the final table
    outfile = paste0(args$out_prefix, '.txt')
    fwrite(table, outfile, sep = '\t')
}

args <- list(
    out_prefix = 'data/genes/220310_common_plofs_to_exclude',
    in_file = 'data/knockouts/alt/ukb_eur_wes_200k_chrCHR_maf0to5e-2_pLoF_damaging_missense.tsv.gz',
    threshold = 10e+4
)

main(args)






