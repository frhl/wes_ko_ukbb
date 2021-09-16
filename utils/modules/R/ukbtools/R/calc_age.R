#' Calculate age
#'
#' By default, calculates the typical "age in years", with a
#' \code{floor} applied so that you are, e.g., 5 years old from
#' 5th birthday through the day before your 6th birthday. Set
#' \code{floor = FALSE} to return decimal ages, and change \code{units}
#' for units other than years.
#' @param time date-of-birth, the day to start calculating age.
#' @param age_day the date on which age is to be calculated.
#' @param units unit to measure age in. Defaults to \code{"years"}. Passed to \link{\code{duration}}.
#' @param floor boolean for whether or not to floor the result. Defaults to \code{TRUE}.
#' @return Age in \code{units}. Will be an integer if \code{floor = TRUE}.
#' @note taken from here \url{https://stackoverflow.com/questions/27096485/change-a-column-from-birth-date-to-age-in-r}
calc_time_since <- function(time, age_day = lubridate::today(), units = "years", floor = TRUE) {
  the_age = lubridate::interval(time, age_day) / lubridate::duration(num = 1, units = units)
  if (floor) return(as.integer(floor(the_age)))
  return(the_age)
}

