#' Simulate a Single Run of a Prey-Predator Model
#'
#' This function simulates a single run of a prey-predator model over a specified
#' period of time, given initial conditions and model parameters.
#'
#' @param n_prey_initial An integer representing the initial number of prey.
#' @param time_max A numeric value indicating the maximum time for the simulation.
#' @param model A function that calculates the propensity given the current number of prey and parameters.
#' @param parameters A list of named parameters required by the model function.
#'
#' @return The number of prey remaining at the end of the simulation.
simulate_single <- function(n_prey_initial, time_max, model, parameters) {
  t <- 0
  n_prey_remaining <- n_prey_initial
  while(t < time_max & n_prey_remaining > 0) {
    r <- stats::runif(1)
    propensity <- model(n_prey_remaining, parameters)

    # Ensure propensity is positive to avoid division by zero or negative time increments
    if (propensity <= 0) {
      warning(paste0("Propensity is non-positive. Parameters are ", parameters))
      break
    }

    tau <- -log(r) / propensity
    t <- t + tau
    n_prey_remaining <- n_prey_remaining - 1
  }

  if(t <= time_max)
    n_prey_remaining
  else
    n_prey_remaining + 1 # at t = time_max, n_prey_remaining was 1 greater
}


#' Simulate Multiple Runs of a Prey-Predator Model
#'
#' This function simulates multiple runs of a prey-predator model, each with the same
#' initial conditions and parameters, and returns a summary of the results.
#'
#' @param n_replicates An integer representing the number of replicate simulations to perform.
#' @inheritParams simulate_single
#'
#' @return A tibble with columns \code{n_prey_initial}, \code{n_prey_eaten}, and \code{n_prey_remaining},
#'         summarizing the results of the simulations. Each row represents a single simulation replicate.
#'
#' @examples
#' # Define a simple linear model for the propensity
#' model <- function(n_prey, params) {
#'   rate <- params$rate
#'   return(rate * n_prey)
#' }
#'
#' # Set the number of replicates, initial number of prey, maximum time, and model parameters
#' n_replicates <- 10
#' n_prey_initial <- 100
#' time_max <- 10
#' parameters <- list(rate = 0.1)
#'
#' # Run the simulations
#' simulate(n_replicates, n_prey_initial, time_max, model, parameters)
#'
#' @export
simulate <- function(n_replicates, n_prey_initial, time_max, model, parameters) {

  if (!is.function(model)) {
    stop("model must be a function")
  }

  if (!is.numeric(n_replicates) || n_replicates <= 0 || round(n_replicates) != n_replicates) {
    stop("n_replicates must be a positive integer")
  }

  if (!is.numeric(n_prey_initial) || n_prey_initial <= 0 || round(n_prey_initial) != n_prey_initial) {
    stop("n_prey_initial must be a positive integer")
  }

  if (!is.numeric(time_max) || time_max <= 0) {
    stop("time_max must be a positive number")
  }

  n_prey_remaining <- vector(length = n_replicates)
  for(i in seq_along(n_prey_remaining)) {
    n_prey_remaining[i] <- simulate_single(n_prey_initial, time_max, model, parameters)
  }

  result <- dplyr::tibble(
    replicate_id = seq_len(n_replicates),
    n_prey_initial = n_prey_initial,
    n_prey_eaten = n_prey_initial - n_prey_remaining,
    n_prey_remaining = n_prey_remaining
  )

  result
}

#' Simulate Multiple Runs of a Prey-Predator Model
#'
#' This function simulates multiple runs of a prey-predator model, each with the same
#' initial conditions and parameters, and returns a summary of the results.
#'
#' @param n_replicates An integer representing the number of replicate simulations to perform.
#' @inheritParams simulate_single
#'
#' @return A tibble with columns \code{n_prey_initial}, \code{n_prey_eaten}, and \code{n_prey_remaining},
#'         summarizing the results of the simulations. Each row represents a single simulation replicate.
#'
#' @examples
#' # Define a simple linear model for the propensity
#' model <- function(n_prey, params) {
#'   rate <- params$rate
#'   return(rate * n_prey)
#' }
#'
#' # Set the number of replicates, initial number of prey, maximum time, and model parameters
#' n_replicates <- 10
#' n_prey_initial <- 100
#' time_max <- 10
#' parameters <- list(rate = 0.1)
#'
#' # Run the simulations
#' simulate_efficient(n_replicates, n_prey_initial, time_max, model, parameters)
#'
#' @export
simulate_efficient <- function(n_replicates, n_prey_initial, time_max, model, parameters) {

  if (!is.function(model)) {
    stop("model must be a function")
  }

  if (!is.numeric(n_replicates) || n_replicates <= 0 || round(n_replicates) != n_replicates) {
    stop("n_replicates must be a positive integer")
  }

  if (!is.numeric(n_prey_initial) || n_prey_initial <= 0 || round(n_prey_initial) != n_prey_initial) {
    stop("n_prey_initial must be a positive integer")
  }

  if (!is.numeric(time_max) || time_max <= 0) {
    stop("time_max must be a positive number")
  }

  # precalculate all possible propensities
  possible_prey_remaining <- n_prey_initial:1
  propensities <- model(possible_prey_remaining, parameters)
  propensities <- rep(propensities, n_replicates)

  # precalculate r ~ unif(0, 1)
  r <- runif(n_replicates * n_prey_initial)

  # calculate taus
  tau <- -log(r) / propensities

  # Reshape tau for easier handling
  tau_matrix <- matrix(tau, nrow = n_prey_initial, ncol = n_replicates, byrow = FALSE)

  # Calculate cumulative sums once for all replicates
  cumsum_matrix <- apply(tau_matrix, 2, cumsum)

  # Find the number of prey remaining
  prey_remaining <- apply(cumsum_matrix, 2, function(times) {
    idx <- which(times > 1)[1]
    if (is.na(idx)) 0 else possible_prey_remaining[idx]
  })

  prey_eaten <- n_prey_initial - prey_remaining
  result <- dplyr::tibble(
    replicate_id = seq_len(n_replicates),
    n_prey_initial = n_prey_initial,
    n_prey_eaten = prey_eaten,
    n_prey_remaining = prey_remaining
  )

  result
}

#' Simulate a Study with Multiple Prey Initial Conditions
#'
#' This function performs simulations for a study based on varying initial prey counts.
#' It uses a provided model and parameters to generate simulation results and combines
#' these results into a single dataframe.
#'
#' @param data A dataframe containing the study parameters. It must include the following columns:
#'   \itemize{
#'     \item \code{prey_initial}: An integer indicating the initial density of prey.
#'     \item \code{n_replicates}: An integer specifying the number of replicates for each \code{prey_initial}.
#'   }
#' @param time_max A positive numeric value specifying the maximum time for the simulation.
#' @param model A function representing the model to be used in the simulation.
#' @param parameters A list of parameters to be passed to the model function.
#'
#' @return A dataframe containing the simulation results.
#'
#' @details
#' The function first validates the inputs to ensure they meet the required specifications.
#' It then iterates over unique values of initial prey density, performs simulations for
#' each value, and combines all results into a single dataframe.
#'
#' If the dataframe has more than one row per unique \code{prey_initial} value or
#' if any other input validations fail, the function will stop with an error message.
#'
#' @examples
#' # Example usage (assuming appropriate model and parameters):
#' data <- data.frame(prey_initial = c(10, 20, 30), n_replicates = c(100, 100, 100))
#' time_max <- 10
#' model <- model_constant()
#' parameters <- list(rate = 0.1)
#' result <- simulate_study(data, time_max, model, parameters)
#' @export
simulate_study <- function(data, time_max, model, parameters) {

  if (!is.data.frame(data) || !("n_prey_initial" %in% colnames(data)) || !("n_replicates" %in% colnames(data))) {
    stop("`data` must be a dataframe containing columns 'n_prey_initial' and 'n_replicates'.")
  }

  if (!is.numeric(time_max) || time_max <= 0) {
    stop("time_max must be a positive number")
  }

  if (!is.function(model)) {
    stop("model must be a function")
  }

  if (!is.list(parameters)) {
    stop("`parameters` must be a list.")
  }


  unique_prey_initial <- sort(unique(data$n_prey_initial))
  if(length(unique_prey_initial) != nrow(data))
    stop("Dataframe should have one row per n_prey_initial.")

  for(i in seq_along(unique_prey_initial)) {

    n_prey_initial <- unique_prey_initial[i]
    tmp <- data %>%
      dplyr::filter(.data$n_prey_initial==n_prey_initial)
    n_replicates <- tmp$n_replicates[1]

    df <- simulate(n_replicates, n_prey_initial, time_max, model, parameters)

    if(i == 1)
      df_stacked <- df
    else
      df_stacked <- df_stacked %>%
      dplyr::bind_rows(df)
  }

  df_stacked
}
