
snp_asGeneticPosLocal <- function (infos.chr, infos.pos, mapdir, ncores = 1, 
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


qc_sumstat <- function(info_snp) {
  # QC summary stats by comparing
  chi2 <- with(info_snp, (beta / beta_se)^2)
  print(round(mean(chi2, na.rm = TRUE), 2))
  S <- rep(NA, ncol(G)); S[info_snp$`_NUM_ID_`] <- chi2
  signif <- pchisq(S, df = 1, lower.tail = FALSE) < 5e-8
  print(sum(signif, na.rm = TRUE))
  ind.keep <- snp_clumping(
    G, infos.chr = map$chr, infos.pos = map$pos, S = S,
    ind.row = ind.val, thr.r2 = 0.01, size = 10e3, ncores = NCORES,
    exclude = which(is.na(S) | !signif))
  print(length(ind.keep))
}

read_hail_sumstat <- function(path_bed){
  # Read summary statistics and format to reference
  # from hail to ldpred2 matching coluns.

  sumstats <- bigreadr::fread2(path_bed)
  # make artificial rsid
  sumstats$rsid <-
      apply(
          data.frame(
              locus = sumstats$locus,
              a0 = sumstats$a0,
              a1 = sumstats$a1
              ),
          1,
          paste,
          collapse = '_')

  # extract position and chromosome
  locus <- do.call(rbind, strsplit(sumstats$locus, split = "\\:"))
  colnames(locus) <- c('chr','pos')
  sumstats <- cbind(sumstats, locus)

  # calculate MAF
  sumstats$MAF <- unlist(lapply(sumstats$AF, function(af) min(af, 1-af)))

  sumstats <-
      data.frame(
          chr = sumstats$chr,
          pos = as.integer(sumstats$pos),
          rsid = sumstats$rsid,
          a1 = sumstats$a1,
          a0 = sumstats$a0,
          n_eff = sumstats$n,
          beta_se = sumstats$standard_error,
          p = sumstats$p_value,
          beta = sumstats$beta,
          INFO = NA,
          MAF = sumstats$MAF
      )

    return(sumstats)
}

file.exists.ext <- function(fname, ext){
    fname <- tools::file_path_sans_ext(fname)
    fname <- paste0(fname, ext)
    return(file.exists(fname))
}



