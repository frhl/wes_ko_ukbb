#' @title SNP as genetic position
#' @param infos.chr vector of chromosomes of length n
#' @param infos.pos vector of genetic positions of length n
#' @param mapdir What directory should genetic maps (interpolated dist) be found at?
#' @param ncores how many cores should be used during computation
#' @return vector of genetic positions of length n
#' @note this is a "hacked" version of florian's function that
#' can work on the cluster.
#' @export


snp_as_genetic_position <- function (infos.chr, infos.pos, mapdir, ncores = 1,
    # get SNP as genetic location (cM) based on either hapmap or omno interpolated distances.
    rsid = NULL, genetic_map = 'omni')
{
    bigassertr::assert_package("R.utils")
    bigassertr::assert_lengths(infos.chr, infos.pos)
    if (!is.null(rsid))
        bigassertr::assert_lengths(rsid, infos.pos)
    snp_split(infos.chr, function(ind.chr, pos, dir, rsid) {
        chr <- attr(ind.chr, "chr")
        stopifnot(genetic_map %in% c('hapmap','omni'))
        map_origin <- ifelse(
            genetic_map %in% 'omni',
            ".OMNI.interpolated_genetic_map",
            ".interpolated_genetic_map")
        basename <- paste0("chr", chr, map_origin)
        gzfile <- file.path(mapdir, paste0(basename, '.gz'))
        mapfile <- file.path(mapdir, basename)
        if (!file.exists(mapfile)) {
            stopifnot(file.exists(gzfile))
            R.utils::gunzip(gzfile)
        }
        map.chr <- bigreadr::fread2(mapfile, showProgress = FALSE)
        if (is.null(rsid)) {
            ind <- bigutilsr::knn_parallel(as.matrix(map.chr$V2),
                as.matrix(pos[ind.chr]), k = 1, ncores = 1)$nn.idx
            new_pos <- map.chr$V3[ind]
        }
        else {
            ind <- match(rsid[ind.chr], map.chr$V1)
            new_pos <- map.chr$V3[ind]
            indNA <- which(is.na(ind))
            if (length(indNA) > 0) {
                pos.chr <- pos[ind.chr]
                new_pos[indNA] <- suppressWarnings(stats::spline(pos.chr,
                  new_pos, xout = pos.chr[indNA], method = "hyman")$y)
            }
        }
        new_pos
    }, combine = "c", pos = infos.pos, dir = dir, rsid = rsid,
        ncores = ncores)
}

