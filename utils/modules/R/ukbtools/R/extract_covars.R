#' @title extract UKBB covaraties
#' @description a function that extracts basisc UKBB parameters
#' @param dt_pheno data.table for phenotypes
#' @param dt_sample_qc data.table for sample QC
#' @param keep what columns should be kept from dt_pheno file?
#' @export

extract_covars <- function(dt_pheno, dt_sample_qc, keep = 'eid'){

  # combine dt and dt_sample_qc
  #keep = c('eid','Submitted.Gender','Inferred.Gender','in.white.British.ancestry.subset')
  remove <- colnames(dt_pheno)[!colnames(dt_pheno) %in% keep]
  dt_sample_qc$eid <- dt_sample_fam$V1
  dt_sample_qc <- dt_sample_qc[,c('eid','Submitted.Gender','Inferred.Gender','in.white.British.ancestry.subset')]
  sum(dt_pheno$eid %in% dt_sample_qc$eid)
  dt <- merge(dt_pheno, dt_sample_qc)

  # get data of birth
  dt$yob <- dt$X34.0.0
  dt$mob <- dt$X52.0.0
  dt$dob <- as.Date(paste0(dt$yob,'-',dt$mob,'-',1))

  # get date of attending recruitment centre
  dt$centre_attend0 <- as.Date(dt$X53.0.0)
  dt$centre_attend1 <- as.Date(dt$X53.1.0)
  dt$centre_attend2 <- as.Date(dt$X53.2.0)

  # get age when attending recruitment centre
  dt$centre_attend_age0 <- calc_time_since(dt$dob, dt$centre_attend0)
  dt$centre_attend_age1 <- calc_time_since(dt$dob, dt$centre_attend1)
  dt$centre_attend_age2 <- calc_time_since(dt$dob, dt$centre_attend2)

  # get age to be used as covariates
  dt$age <- dt$centre_attend_age0
  dt$age1 <- dt$centre_attend_age0^2
  dt$age2 <- dt$centre_attend_age0^3

  # get principal components
  bool_pc <- grep(22009,colnames(dt))
  n <- length(colnames(dt)[bool_pc])
  colnames(dt)[bool_pc] <- paste0('PC',1:n)

  # negative coding from BiLEVE (recoded 1); positive coding from axiom (recoded 2)
  dt$genotyping.array <- ifelse(dt$X22000.0.0 <0, 1, 2)
  dt$sex = ifelse(dt$Inferred.Gender == "F", 2, ifelse(dt$Inferred.Gender == "M", 1, NA))
  colnames(dt)[colnames(dt) == "X54.0.0"] = "ukbb.centre"

  # clean up
  bool_keep <- ! colnames(dt) %in% remove
  dt = dt[,bool_keep, with = FALSE]

  return(dt)
}
