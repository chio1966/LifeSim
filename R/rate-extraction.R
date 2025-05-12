#' Calculate Age in Years
#'
#' Computes the age (as a fractional number of years) between two dates.
#'
#' @param birthdate A `Date` object representing the date of birth.
#' @param age_date A `Date` object representing the later date at which age is to be calculated.
#'
#' @return A numeric value giving the age in years (including decimal fraction).
#'
#' @examples
#' library(lubridate)
#' birthdate <- ymd("1980-06-15")
#' age_date <- ymd("2020-06-14")
#' age_function(birthdate, age_date)
age_function <- function(birthdate,age_date) {
  return(time_length(interval(birthdate, age_date), "years"))
}

#' Adjust Leap Day (Feb 29) to Feb 28
#'
#' This function checks for any dates that fall on February 29 and adjusts them back to February 28,
#' to avoid issues when performing date calculations (e.g., in non-leap years).
#'
#' @param date A `Date` vector.
#'
#' @return A `Date` vector where all leap days (Feb 29) are moved to Feb 28.
#'
adjust_elipe <- function(date) {
  # Move Feb 29 to Feb 28
  date[month(date) == 2 & day(date) == 29] <- date[month(date) == 2 & day(date) == 29] - days(1)

  return(date)
}

#' Extract Life Table Rates and Time Intervals
#'
#' Given a data frame of individuals with diagnosis dates, demographic covariates, and a life table,
#' this function computes a list of time intervals and corresponding hazard rates for each individual.
#' Useful for modeling survival using population-level expected mortality.
#'
#' @param df A `data.frame` containing individual-level data. Must include columns for age, sex, etc., as described in `mapping`.
#' @param ratetable A multi-dimensional life table object, such as those from `survival::survexp.us`.
#' @param mapping A named character vector mapping column names in `df` to dimension names in `ratetable`. The first two elements must be the names corresponding to age and year.
#'   E.g., `c("age" = "age", "year" = "year", "sex" = "sex", "ethnicity" = "race")`.
#' @param start A string giving the column name in `df` for diagnosis or start date.
#' @param end A string giving the column name in `df` for end of follow-up.
#' @param rate_scale A numeric multiplier to rescale extracted rates this will depends on the ratetable.
#' @param time_scale A numeric multiplier to rescale time units the default time is given in days, adjust accordingly.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{times}{A list of numeric vectors giving time intervals per individual in days}
#'   \item{rates}{A list of numeric vectors giving hazard rates per interval per individual.}
#' }
#' @examples
#' library(lubridate)
#' library(dplyr)
#' library(survival)
#'
#' # Simulate test data
#' set.seed(123)
#' n <- 200
#'
#' # Random diagnosis dates between 1970 and 1972
#' start_date <- as.Date("1970-01-01") + runif(n, 0, as.numeric(as.Date("1972-12-31") - as.Date("1970-01-01")))
#'
#' # Random end dates between 1985 and 1986
#' end_date <- as.Date("1985-01-01") + runif(n, 0, as.numeric(as.Date("1986-12-31") - as.Date("1985-01-01")))
#'
#' # Ages at diagnosis between 40 and 70
#' age <- runif(n, 40, 70)
#'
#' # Sex and ethnicity as factors
#' sex <- sample(c("male", "female"), n, replace = TRUE)
#' ethnicity <- sample(c("white", "black"), n, replace = TRUE)
#'
#' # Assemble data frame
#' test_df <- data.frame(
#'   start_date = start_date,
#'   end_date = end_date,
#'   age = age,
#'   sex = sex,
#'   ethnicity = ethnicity
#' )
#'
#' # Define mapping to ratetable dimensions
#' mapping <- c("age" = "age", "year" = "year", "sex" = "sex", "ethnicity" = "race")
#'
#' # Use the population rate table from the 'survival' package
#' ratetable <- survival::survexp.usr
#'
#' # Run the extraction
#' result <- rate_extraction(
#'   df = test_df,
#'   ratetable = ratetable,
#'   mapping = mapping,
#'   start = "start_date",
#'   end = "end_date",
#'   rate_scale = 365.25,      # Convert yearly rates to daily
#'   time_scale = 1 / 365.25   # Convert days to years
#' )
#'
#' # Inspect structure
#' str(result)
#'
#'
rate_extraction <- function(df, ratetable, mapping, start, end, rate_scale, time_scale) {

times <- list()
rates <- list()

for ( i in seq_len(nrow(df))){

  df_temp <- df[i, ]

  #start and end will be supplied as arguments
  diagnosis_date <- adjust_elipe(df_temp[[start]])
  end_date <- adjust_elipe(df_temp[[end]])
  age_at_diagnosis <- df_temp[[names(mapping)[1]]]

  dob <- diagnosis_date - years(floor(age_at_diagnosis))
  dob <- dob - days(floor((age_at_diagnosis-floor(age_at_diagnosis))*365.25))

  dob <- adjust_elipe(dob)

  jan_1st_dates <- seq(from = floor_date(diagnosis_date, "year") + years(1),
                       to = end_date, by = "year")

  jan_1st_dates <- jan_1st_dates[jan_1st_dates >= diagnosis_date]

  birthday <- seq(from = make_date(year(diagnosis_date), month(dob), day(dob)), to = end_date, by = "year")
  birthday <- birthday[birthday >= diagnosis_date]

  dates <- unique(c(diagnosis_date, end_date, jan_1st_dates, birthday))
  dates <- sort(dates)
  #now need to make function for ages at those dates ect

  age_at_dates <- floor(age_function(birthdate = dob, age_date = dates))
  year_at_dates <- year(dates)

  factors <- df_temp %>% select(names(mapping)[3:length(mapping)]) %>% slice(rep(1,length(age_at_dates)))

  data_temp <- bind_cols(data.frame(age_at_dates,year_at_dates),factors)

  names(data_temp)[1:2] <- names(mapping)[1:2] #so now data_temp has names = names(mapping)



  columns_order <- match(names(dimnames(ratetable)),mapping[names(data_temp)]) #reordering df columns to match the dims of ratetable

  data_temp_ordered <- data_temp[,columns_order]
  ##
  rates_temp <- numeric(nrow(data_temp_ordered))

  for (j in 1:nrow(data_temp_ordered)) {
    rate_table_arg <- as.list(as.character(data_temp_ordered[j, ]))
    rates_temp[j] <- do.call(`[`, c(list(ratetable), rate_table_arg))
  }
  times[[i]] <- c(0,diff(dates)) * time_scale
  rates[[i]] <- rates_temp * rate_scale
}
return(list(times = times, rates = rates))
}
