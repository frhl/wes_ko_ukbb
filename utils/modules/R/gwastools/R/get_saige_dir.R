#' @title list directory for SAIGE analysis
#' @param cond either "none", "common", "rare", "combined"
#' @return a string for saige dir

get_saige_dir <- function(directory = "") {
    if (directory == "") {
        out_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2/min_mac4"
    }
    else if (directory == "common") {
        out_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_common/min_mac4"
    }
    else if (directory == "rare") {
        out_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_rare_cond/min_mac4"
    }
    else if (directory == "combined") {
        out_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_collapsed_urv/min_mac4"
    }
    else if (directory == "chet_only") {
        out_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_only_chets/min_mac4"
    }
    else if (directory == "hom_only") {
        out_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_only_homs/min_mac4"
    }
    else if (directory == "add_encoding") {
        out_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_encoding_012/min_mac4"
    }
    else if (directory == "only_singletons") {
        out_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_only_singletons/min_mac4"
    }
    else if (directory == "exclude_singletons") {
        out_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_exclude_singletons/min_mac4"
    }
    else {
        stop(paste(directory, "is not valid saige directory subset"))
    }
    if (!dir.exists(out_dir)){
        stop(paste(out_dir, "is not a directory!"))
    }
    return(out_dir)
}


