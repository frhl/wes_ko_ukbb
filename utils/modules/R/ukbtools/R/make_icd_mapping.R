#' @title make icd mapping funcs
#' @description Create a function that allows quick mapping between ICD and disease cateogories
#' @param from a string
#' @param to a string
#' @param icd_path path to icd bridge file
#' @export
#' @examples
#' \dontrun{
#' f <- create_icd_mapping('ICD10','disease')
#' f(c('124','Z851'))
#' }

make_icd_mapping <- function(from, to, icd_path = "/well/lindgren/UKBIOBANK/ferreira/convert_disease_codes/icd9_icd10_lkps.txt"){

  # Read icd code
  icd <- fread(icd_path)
  colnames(icd) <- gsub('_code','',colnames(icd))
  stopifnot(c(from,to) %in% colnames(icd))

  # Deal with NAs
  icd <- icd[,colnames(icd) %in% c(from, to),with = F]
  icd <- icd[!duplicated(icd)]
  icd[icd==''] <- NA
  icd[icd=='UNDEF'] <- NA
  icd <- icd[!is.na(icd[[from]]), ]
  mapping <- as.list(icd[[to]])
  names(mapping) <- icd[[from]]

  # return function that can do the mapping
  mapper <- function(code){
    map <- mapping[code]
    ids <- names(map)
    is_na <- is.na(ids)
    names(map) <- code
    map[is_na] <- NA
    return(as.character(unlist(map)))
  }

  return(mapper)
}




