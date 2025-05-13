#' Simulate survival times based on life table rates
#'
#' This function simulates survival times and censoring statuses based on a given life table
#' with cumulative hazard values, using a uniform random sample and inverse cumulative hazard function.
#'
#' @param n An integer specifying the number of simulations (individuals).
#' @param lst A list containing two elements: `times` and `rates`.
#' `times` is a list of time intervals, and `rates` is a list of corresponding rates for each individual.
#' @param seed An integer specifying the seed for random number generation for reproducibility.
#'
#' @returns A list containing:
#'   \item{time}{A list of simulated survival times for each individual.}
#'   \item{status}{A list of censoring statuses (1 for event, 0 for censored) for each individual.}
#'
#' @export
#'
#' @examples
#' library(survival)
#' library(dplyr)
#' library(lubridate)
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

lifetable_sim = function(n, lst, seed) {

  if (missing(seed) || !is.numeric(seed) || length(seed) != 1 || is.na(seed)) {
    stop("Argument 'seed' must be a single numeric value and cannot be missing.")
  }

  # Check 'n'
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
    stop("Argument 'n' must be a single positive integer.")
  }

  set.seed(seed)
  time = list()
  status = list()

  for (i in seq_along(lst$times)) {
    ts = cumsum(lst$times[[i]])
    haz = lst$rates[[i]]

    # Cumulative hazard function values at ts.
    # Cumulative hazard function of a piecewise constant hazard function is piecewise linear function.
    H_values = cumsum(diff(ts)* haz[-length(haz)])
    Hinv = approxfun(x = c(0, H_values), y = ts, method = "linear")

    # Simulate uniform sample for inversion
    uniform_sample = -log(runif(n))

    # preallocating storage vectors
    time_temp = numeric(n)
    status_temp = integer(n)
    # Index's of any uniform sample that would constitute a survival time greater than the observation window
    is_censored = uniform_sample > max(H_values)

    # If TRUE then they are censored at the length of the observation window
    time_temp[is_censored] = max(ts)

    # If FALSE then the survival time is obtained by inverse transform sampling method
    time_temp[!is_censored] = Hinv(uniform_sample[!is_censored])

    # Assign status: 0 = censored, 1 = event
    status_temp[is_censored] = 0
    status_temp[!is_censored] = 1

    # Storing results in main list
    time[[i]] = time_temp
    status[[i]] = status_temp
  }
  return(list(time = time, status = status))
}
