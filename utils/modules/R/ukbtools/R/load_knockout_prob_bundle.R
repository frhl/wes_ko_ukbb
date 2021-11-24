load_knockout_prob_bundle <- function(files, maf, mutation){
    
    # subset to files
    bool_maf = grepl(maf, files)
    mutation_file = paste0(mutation,'_','ko_prob')
    bool_mutation = grepl(mutation_file, files)

    # load files
    files <- files[bool_mutation & bool_maf]
    stopifnot(length(files) > 0)
    dt <- setDT(do.call(rbind, lapply(files, zcat)))

    return(dt)
}



