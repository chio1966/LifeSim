#' Create Individual Hazard Function Objects
#'
#' This function generates a list of hazard (`h`), cumulative hazard (`H`), and inverse cumulative hazard (`Hinv`)
#' functions for each individual, based on piecewise constant rates over specified intervals.
#'
#' @param lst A list containing `times` and `rates`, typically the output from `rate_extraction()`.
#'
#' @returns A list where each element is a list of three functions (`h`, `H`, and `Hinv`) for each individual:
#' * `h`: the hazard function over time.
#' * `H`: the cumulative hazard function.
#' * `Hinv`: the inverse cumulative hazard function (for sampling survival times).
#'
#' @export
#'
#' @examples
#' library(survival)
#' library(lubridate)
#' library(dplyr)
#' # Generate example input data
#' set.seed(123)
#' n <- 200
#' start_date <- as.Date("1970-01-01") + runif(n, 0, as.numeric(as.Date("1972-12-31") - as.Date("1970-01-01")))
#' end_date <- as.Date("1985-01-01") + runif(n, 0, as.numeric(as.Date("1986-12-31") - as.Date("1985-01-01")))
#' age <- runif(n, 40, 70)
#' sex <- sample(c("male", "female"), n, replace = TRUE)
#' ethnicity <- sample(c("white", "black"), n, replace = TRUE)
#'
#' test_df <- data.frame(
#'   start_date = start_date,
#'   end_date = end_date,
#'   age = age,
#'   sex = sex,
#'   ethnicity = ethnicity
#' )
#'
#' mapping <- c("age" = "age", "year" = "year", "sex" = "sex", "ethnicity" = "race")
#'
#' # Extract lifetable times and rates
#' rates_out <- rate_extraction(test_df, ratetable = survexp.usr, mapping = mapping,
#'                               start = "start_date", end = "end_date",
#'                               rate_scale = 365.25, time_scale = 1/365.25)
#'
#' # Simulate survival times based on extracted rates
#' result <- lifetable_sim(n = 1, lst = rates_out, seed = 42)
#' str(result)

individual_functions = function(lst){
  function_list <- list()
  for( i in seq_along(lst$times)) {

    ts <- cumsum(lst$times[[i]])
    haz <- lst$rates[[i]]
    H_values <- cumsum(diff(ts)* haz[-length(haz)])

    function_list[[i]] <- list(
      h = approxfun(x = ts, y = haz, method = "constant"),
      H = approxfun(x = ts, y = c(0, H_values), method = "linear"),
      Hinv = approxfun(x = c(0, H_values), y = ts, method = "linear")
    )
  }
  return(function_list)
}
