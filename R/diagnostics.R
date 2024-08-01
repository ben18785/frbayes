
#' Create Bootstrapped Samples
#'
#' This function generates bootstrapped samples from a given dataset by simulating
#' a study multiple times using a specified model and parameters.
#'
#' @param n_bootstraps An integer representing the number of bootstrap samples to generate.
#' @param data A data frame containing the experimental data, which should include a column `n_prey_initial`.
#' @param time_max A numeric value indicating the maximum time for the simulation.
#' @param model A function that calculates the propensity given the current number of prey and parameters.
#' @param mle_parameters A list of maximum likelihood estimated parameters required by the model function.
#'
#' @return A data frame containing the bootstrapped samples. The data frame includes columns
#'   `n_prey_initial`, `n_replicates`, and `bootstrap_id` for each bootstrap sample.
#'
#' @import dplyr
create_bootstrapped_samples <- function(
    n_bootstraps,
    data,
    time_max,
    model,
    mle_parameters) {

  experimental_setup <- data %>%
    group_by(.data$n_prey_initial) %>%
    count() %>%
    dplyr::rename(n_replicates=.data$n)

  for(i in seq(1, n_bootstraps, 1)) {

    df_sim <- simulate_study(
      data=experimental_setup,
      time_max = time_max,
      model = model,
      parameters = mle_parameters
    ) %>%
      dplyr::mutate(
        bootstrap_id=i
      )

    if(i == 1)
      df_bootstrap_samples <- df_sim
    else
      df_bootstrap_samples <- df_bootstrap_samples %>%
        dplyr::bind_rows(df_sim)

  }

  df_bootstrap_samples
}

#' Create Empirical CDFs for Prey Eaten Data
#'
#' This function computes the empirical cumulative distribution functions (ECDFs)
#' for the number of prey eaten, grouped by the initial number of prey.
#'
#' @param data A data frame containing the experimental data, which should include columns
#'   \code{n_prey_initial} and \code{n_prey_eaten}.
#'
#' @return A data frame with columns \code{n_prey_eaten}, \code{ecdf}, and \code{n_prey_initial}.
#'
#' @details
#' The function calculates the ECDF of the number of prey eaten for each unique value
#' of \code{n_prey_initial} in the input data. The ECDF is computed over a range of possible
#' values for the number of prey eaten, from 0 to the initial number of prey.
#' @import dplyr
create_study_ecdfs <- function(data) {

  prey_initial <- sort(unique(data$n_prey_initial))

  for(i in seq_along(prey_initial)) {

    prey_initial_tmp <- prey_initial[i]

    df_short <- data %>%
      dplyr::filter(.data$n_prey_initial==prey_initial_tmp)

    # calculate ecdf series
    f_ecdf <- stats::ecdf(df_short$n_prey_eaten)
    prey_eaten_range <- seq(0, prey_initial_tmp, 1)
    series_ecdf <- f_ecdf(prey_eaten_range)

    df_tmp <- dplyr::tibble(
      n_prey_eaten=prey_eaten_range,
      ecdf=series_ecdf,
      n_prey_initial=prey_initial_tmp
    )

    if(i == 1)
      df_ecdfs <- df_tmp
    else
      df_ecdfs <- df_ecdfs %>% dplyr::bind_rows(df_tmp)
  }

  df_ecdfs
}

#' Create Bootstrapped ECDFs for Simulated Data
#'
#' This function computes empirical cumulative distribution functions (ECDFs)
#' for bootstrapped samples of simulated data.
#'
#' @param df_bootstrap_samples A data frame containing bootstrapped samples with columns
#'   \code{bootstrap_id}, \code{n_prey_initial}, and \code{n_prey_eaten}.
#'
#' @return A data frame with columns \code{n_prey_eaten}, \code{ecdf_sim},
#'   \code{n_prey_initial}, and \code{bootstrap_id}, representing the ECDF values
#'   for each bootstrapped sample.
#'
#' @details
#' The function calculates the ECDF of the number of prey eaten for each unique value
#' of \code{n_prey_initial} in the bootstrapped samples. The ECDF is computed for each
#' bootstrap replicate separately, and the results are combined into a single data frame.
#' @import dplyr
create_bootstrapped_ecdf_simulated <- function(
    df_bootstrap_samples
) {

  n_bootstraps <- max(df_bootstrap_samples$bootstrap_id)
  for(i in 1:n_bootstraps) {

    df_bootstrap_samples_tmp <- df_bootstrap_samples %>%
      dplyr::filter(.data$bootstrap_id==i)

    df_ecdfs <- create_study_ecdfs(df_bootstrap_samples_tmp) %>%
      dplyr::mutate(bootstrap_id=i)

    if(i == 1)
      df_ecdfs_sim <- df_ecdfs
    else
      df_ecdfs_sim <- df_ecdfs_sim %>% dplyr::bind_rows(df_ecdfs)
  }
  df_ecdfs_sim <- df_ecdfs_sim %>%
    dplyr::rename(ecdf_sim=.data$ecdf)

  df_ecdfs_sim
}


#' Create Bootstrapped ECDFs for Real and Simulated Data
#'
#' This function computes empirical cumulative distribution functions (ECDFs)
#' for both real and bootstrapped simulated data, and combines the results into
#' a single data frame.
#'
#' @param n_bootstraps An integer representing the number of bootstrap samples to generate.
#' @param data A data frame containing the real data with columns \code{n_prey_initial} and \code{n_prey_eaten}.
#' @param time_max A numeric value indicating the maximum time for the simulation.
#' @param model A function that calculates the propensity given the current number of prey and parameters.
#' @param mle_parameters A list of named parameters required by the model function.
#'
#' @return A data frame with columns \code{n_prey_eaten}, \code{ecdf_sim},
#'   \code{ecdf_real}, \code{n_prey_initial}, and \code{bootstrap_id}, representing the ECDF values
#'   for both real and bootstrapped simulated data.
#'
#' @details
#' The function first generates bootstrapped samples of the simulated data, computes
#' the ECDF of the number of prey eaten for each unique value of \code{n_prey_initial}
#' in both the real and simulated data, and then combines the results into a single
#' data frame for comparison.
#'
#' @examples
#' # Example real data
#' experimental_setup <- data.frame(n_prey_initial = c(10, 20, 30),
#' n_replicates = c(100, 100, 100))
#' time_max <- 10
#' model <- model_stochastic_degradation()
#' parameters <- list(rate = 0.1)
#' data <- simulate_study(experimental_setup, time_max, model, parameters)
#'
#' # Suppose this maximum likelihood estimate
#' mle_parameters <- list(rate = 0.1)
#'
#' # Create bootstrapped ECDFs for real and simulated data
#' ecdfs_both <- create_bootstrapped_ecdf_real_simulated(
#'   n_bootstraps = 10,
#'   data = data,
#'   time_max = time_max,
#'   model = model,
#'   mle_parameters = mle_parameters
#' )
#'
#' @import dplyr
#' @export
create_bootstrapped_ecdf_real_simulated <- function(
    n_bootstraps,
    data,
    time_max,
    model,
    mle_parameters) {

  # Input validation
  if (!is.numeric(n_bootstraps) || length(n_bootstraps) != 1 || n_bootstraps <= 0 || n_bootstraps != as.integer(n_bootstraps)) {
    stop("Parameter 'n_bootstraps' must be a positive integer.")
  }

  if (!is.data.frame(data)) {
    stop("Parameter 'data' must be a data frame.")
  }

  if (!all(c("n_prey_initial", "n_prey_eaten") %in% names(data))) {
    stop("Data frame must contain columns 'n_prey_initial' and 'n_prey_eaten'.")
  }

  if (!is.numeric(time_max) || length(time_max) != 1 || time_max <= 0) {
    stop("Parameter 'time_max' must be a positive numeric value.")
  }

  if (!is.function(model)) {
    stop("Parameter 'model' must be a function.")
  }

  if (!is.list(mle_parameters) || any(sapply(mle_parameters, is.null))) {
    stop("Parameter 'mle_parameters' must be a non-null list of named parameters.")
  }

  df_bootstrap_samples <- create_bootstrapped_samples(
    n_bootstraps,
    data,
    time_max,
    model,
    mle_parameters
  )

  df_ecdfs_sim <- create_bootstrapped_ecdf_simulated(df_bootstrap_samples)

  df_ecdfs_real <- create_study_ecdfs(data) %>%
    dplyr::rename(ecdf_real=.data$ecdf)

  df_ecdfs_both <- df_ecdfs_sim %>%
    dplyr::left_join(df_ecdfs_real, by=c("n_prey_eaten", "n_prey_initial"))

  df_ecdfs_both
}
