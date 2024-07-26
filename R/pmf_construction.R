

#' Create Approximate Counts from Simulation Results
#'
#' This function processes a dataframe of simulation results to count occurrences of each
#' value of `n_prey_remaining`. It ensures that all possible values of `n_prey_remaining`
#' from 0 to `n_prey_initial` are represented, even if they were not observed in the
#' simulation. It replaces missing values with zero and provides a complete dataframe
#' with counts for each possible value.
#'
#' @param simulation_result A dataframe containing simulation results with a column
#'   `n_prey_remaining` indicating the number of prey remaining.
#' @param n_prey_initial An integer representing the initial number of prey. It defines
#'   the range of possible values for `n_prey_remaining`.
#'
#' @return A tibble with two columns:
#' \item{n_prey_remaining}{An integer indicating the number of prey remaining.}
#' \item{n}{An integer indicating the count of occurrences for each value of
#'   `n_prey_remaining`. Missing values are replaced with zero.}
#'
#' @examples
#' # Example simulation result
#' simulation_result <- tibble(n_prey_remaining = c(0, 1, 1, 2, 2, 2, 3))
#'
#' # Generate approximate counts
#' create_approx_counts(simulation_result, n_prey_initial = 3)
#'
#' @import dplyr
#' @import tidyr
create_approx_counts <- function(simulation_result, n_prey_initial) {

  # Count occurrences of each value of n_prey_remaining
    df_counts <- simulation_result %>%
      dplyr::count(n_prey_remaining)

  # Create a tibble with all possible values of n_prey_remaining
  df_all <- dplyr::tibble(n_prey_remaining = 0:n_prey_initial)

  # Merge counts with the tibble of all possible values
  df_all %>%
    dplyr::left_join(df_counts, by = "n_prey_remaining") %>%
    tidyr::replace_na(list(n = 0)) %>%
    dplyr::arrange(n_prey_remaining)
}

#' Construct Empirical PMF with Log Probabilities
#'
#' This function constructs a dataframe containing the empirical probability mass function (PMF)
#' and its logarithm from simulation results. It computes the PMF using a Dirichlet Process
#' posterior mean, with a uniform prior over all possible values of `n_prey_remaining`. The
#' resulting dataframe includes the PMF and its logarithm for each possible value of `n_prey_remaining`.
#'
#' @inheritParams create_approx_counts
#' @param alpha A numeric value representing the hyperparameter for the Dirichlet Process prior.
#'   It controls the influence of the prior in the posterior calculation.
#'
#' @return A tibble with the following columns:
#' \item{n_prey_remaining}{An integer indicating the number of prey remaining.}
#' \item{pmf}{The empirical probability mass function for each value of `n_prey_remaining`.}
#' \item{log_prob}{The logarithm of the PMF.}
#'
#' @details
#' The PMF is computed using the Dirichlet Process posterior mean, where the prior measure is
#' uniform across all possible values of `n_prey_remaining`. The PMF is adjusted by the hyperparameter
#' `alpha`, and its logarithm is calculated for each value.
#'
#' @examples
#' # Example usage with hypothetical simulation results and parameters
#' simulation_result <- tibble(n_prey_remaining = c(0, 1, 1, 2, 2, 2, 3))
#' n_prey_initial <- 3
#' alpha <- 1
#'
#' # Construct the empirical PMF with log probabilities
#' df_pmf_log <- construct_empirical_pmf_log_df(simulation_result, n_prey_initial, alpha)
#' print(df_pmf_log)
#'
#' @import dplyr
#' @importFrom tidyr replace_na
construct_empirical_pmf_log_df <- function(simulation_result, n_prey_initial, alpha) {

  # prior measure is uniform on set of possible n_prey_remaining
  p <- 1 / (n_prey_initial + 1)
  H <- rep(p, (n_prey_initial + 1))

  # use Dirichlet Process posterior mean as posterior measure
  df_pmf <- create_approx_counts(simulation_result, n_prey_initial) %>%
    dplyr::mutate(H=H) %>%
    dplyr::mutate(pmf=alpha*H+n) %>%
    dplyr::mutate(pmf=pmf/sum(pmf)) %>%
    dplyr::select(-c("n", "H")) %>%
    dplyr::mutate(log_prob=log(pmf))

  df_pmf
}

#' Create a Log Probability Mass Function (PMF)
#'
#' Constructs a function that computes the logarithm of the probability mass function (PMF)
#' for a given number of prey remaining based on empirical data from simulations and a Dirichlet
#' process prior.
#'
#' This function generates a log probability mass function (PMF) by first calculating the
#' empirical PMF from simulation results and then using it to compute the logarithm of the PMF
#' for specified numbers of prey remaining. The function handles valid and invalid inputs and
#' provides appropriate outputs.
#'
#' @inheritParams construct_empirical_pmf_log_df
#'
#' @return A function that takes a single argument, `n_prey_remaining`, and returns the logarithm
#'   of the PMF for that value. If the input is out of bounds or not a non-negative integer, the
#'   function returns `-Inf`.
#'
#' @examples
#' # Example data frame from simulations
#' simulation_result <- tibble::tibble(n_prey_remaining = c(0, 1, 2, 1, 0, 2))
#' n_prey_initial <- 2
#' alpha <- 1
#'
#' # Create the log PMF function
#' pmf_log_function <- create_pmf_log(simulation_result, n_prey_initial, alpha)
#'
#' # Use the log PMF function
#' log_prob_0 <- pmf_log_function(0)
#' log_prob_1 <- pmf_log_function(1)
#' log_prob_2 <- pmf_log_function(2)
#'
#' # Check behavior for out-of-bounds values
#' log_prob_neg <- pmf_log_function(-1)
#' log_prob_out_of_bound <- pmf_log_function(3)
create_pmf_log <- function(simulation_result, n_prey_initial, alpha) {

  df_pmf <- construct_empirical_pmf_log_df(simulation_result, n_prey_initial, alpha)

  pmf_log <- function(n_prey_remaining) {

    if(n_prey_remaining < 0 || n_prey_remaining > n_prey_initial)
      return(-Inf)

    if (!is.numeric(n_prey_remaining) || round(n_prey_remaining) != n_prey_remaining) {
      return(-Inf)
    }

    df_pmf$log_prob[n_prey_remaining + 1]
  }

  pmf_log
}

#' Calculate Log Probability of Multiple Prey Remaining States
#'
#' Computes the log probability of observing specific numbers of prey remaining in a simulation
#' based on a model, given an initial number of prey and other parameters. This function uses
#' simulated data to estimate the log probability of each state using a probability mass function
#' (PMF) and aggregates these probabilities for a list of prey remaining states.
#'
#' The function first generates simulated data using a given model and parameters, then creates
#' a log PMF function from this data. It calculates the log probability for each specified number
#' of prey remaining by summing the log PMF values.
#'
#' @param parameters A list or vector containing parameters required by the model function.
#'   The specific parameters depend on the model used.
#' @param ns_prey_remaining A numeric vector of integers representing the different states of
#'   prey remaining for which the log probability needs to be calculated.
#' @param n_prey_initial An integer specifying the initial number of prey. This is used to define
#'   the range of possible values for prey remaining.
#' @param model A function that computes the propensity or rate of change based on the number of
#'   prey remaining and the parameters provided. The model should be compatible with the simulation
#'   function.
#' @param time_max A numeric value representing the maximum time for the simulation.
#' @param n_replicates An integer specifying the number of simulation runs to perform. Default is 1000.
#' @param alpha A numeric parameter for the Dirichlet process prior. Default is 1.
#'
#' @return A numeric value representing the total log probability of observing the specified
#'   states of prey remaining. This is the sum of the log probabilities for each value in
#'   `ns_prey_remaining`.
#'
#' @examples
#' # Example usage
#' parameters <- list(rate = 0.1)  # Example parameter for a constant rate model
#' ns_prey_remaining <- c(0, 1, 2)  # States of prey remaining
#' n_prey_initial <- 2
#' model <- model_constant_rate()  # Using a constant rate model
#' time_max <- 10
#'
#' # Compute log probabilities
#' log_prob <- log_probability_single_prey_initial(
#'   parameters = parameters,
#'   ns_prey_remaining = ns_prey_remaining,
#'   n_prey_initial = n_prey_initial,
#'   model = model,
#'   time_max = time_max,
#'   n_replicates = 1000,
#'   alpha = 1
#' )
#'
log_probability_single_prey_initial <- function(
    parameters,
    ns_prey_remaining, n_prey_initial,
    model, time_max,
    n_replicates=1000, alpha=1) {

  simulation_result <- simulate(n_replicates, n_prey_initial, time_max, model, parameters)
  log_prob_function <- create_pmf_log(simulation_result, n_prey_initial, alpha)

  log_prob <- 0
  for(i in seq_along(ns_prey_remaining)) {
    log_prob = log_prob + log_prob_function(ns_prey_remaining[i])
  }

  log_prob
}


log_probability <- function(
    parameters,
    data,
    model,
    time_max=1, n_replicates=1000, alpha=1) {

  # Validate inputs
  if (!is.data.frame(data) && !is_tibble(data)) {
    stop("`data` must be a dataframe or tibble.")
  }

  required_columns <- c("n_prey_initial", "n_prey_remaining")
  missing_columns <- setdiff(required_columns, colnames(data))

  if (length(missing_columns) > 0) {
    stop("`data` must contain the following columns: ", paste(missing_columns, collapse = ", "))
  }

  if (!is.numeric(time_max) || time_max <= 0) {
    stop("`time_max` must be a positive number.")
  }

  if (!is.numeric(n_replicates) || n_replicates <= 0 || round(n_replicates) != n_replicates) {
    stop("`n_replicates` must be a positive integer.")
  }

  if (!is.numeric(alpha) || alpha <= 0) {
    stop("`alpha` must be a positive number.")
  }

  log_prob <- 0
  unique_prey_initial <- sort(unique(data$n_prey_initial))
  for(i in seq_along(unique_prey_initial)) {

    df_single_prey_initial <- data %>%
      dplyr::filter(n_prey_initial==unique_prey_initial[i])

    log_prob_increment <- log_probability_single_prey_initial(
      parameters = parameters,
      ns_prey_remaining = df_single_prey_initial$n_prey_remaining,
      n_prey_initial = unique_prey_initial[i],
      model = model,
      time_max = time_max,
      n_replicates = n_replicates,
      alpha = alpha
    )

    log_prob <- log_prob + log_prob_increment

  }

  log_prob
}
