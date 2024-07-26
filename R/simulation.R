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
    r <- runif(1)
    propensity <- model(n_prey_remaining, parameters)

    # Ensure propensity is positive to avoid division by zero or negative time increments
    if (propensity <= 0) {
      warning("Propensity is non-positive. Stopping simulation.")
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

#' Generate a Constant Rate Model for Prey-Predator Simulation
#'
#' This function returns a model function where the rate of prey consumption
#' is constant and proportional to the current number of prey. The returned function
#' can be used in simulations of prey-predator interactions.
#'
#' @return A function that computes the propensity of prey consumption based on the current
#' number of prey and a set of parameters. The returned function takes two arguments:
#' \describe{
#'   \item{prey}{Integer representing the current number of prey.}
#'   \item{parameters}{A list containing the parameters for the model. Must include:}
#'     \describe{
#'       \item{rate}{A numeric value representing the constant rate of prey consumption per prey.}
#'     }
#' }
#'
#' @examples
#' model <- model_constant_rate()
#' parameters <- list(rate = 0.1)
#' model(10, parameters)  # Returns 1.0
#'
#' @export
model_constant_rate <- function() {
  function(prey, parameters) {
    rate <- parameters$rate
    rate * prey
  }
}


# from supplementary here: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13039

#' Generate a Rogers II Model Function for Prey-Predator Simulation
#'
#' This function returns a model function based on the Rogers II model, where the rate of prey consumption
#' is described by a functional response that depends on the current number of prey and two parameters:
#' `a` (the maximum rate of consumption) and `h` (the handling time).
#'
#' The Rogers II model is defined as:
#' \deqn{ \text{rate} = \frac{a \cdot \text{prey}}{1 + a \cdot h \cdot \text{prey}} }
#' where \code{a} is the maximum rate of consumption and \code{h} is the handling time.
#'
#' @return A function that calculates the propensity of prey consumption based on the current
#' number of prey and a set of parameters. The returned function takes two arguments:
#' \describe{
#'   \item{prey}{Numeric representing the current number of prey.}
#'   \item{parameters}{A list containing the parameters for the model. Must include:}
#'     \describe{
#'       \item{a}{A numeric value representing the maximum rate of consumption.}
#'       \item{h}{A numeric value representing the handling time.}
#'     }
#' }
#'
#' @examples
#' model <- model_rogersII()
#' parameters <- list(a = 0.2, h = 0.5)
#' model(10, parameters)  # Computes the rate of prey consumption for 10 prey
#'
#' @export
model_rogersII <- function() {
  function(prey, parameters) {
    a <- parameters$a
    h <- parameters$h
    a * prey / (1 + a * h * prey)
  }
}

# from supplementary here: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13039

#' Generate a Type III Model Function for Prey-Predator Simulation
#'
#' This function returns a model function based on the Type III functional response model,
#' where the rate of prey consumption is described by a quadratic function of the
#' number of prey and two parameters: `b` (the maximum rate of consumption) and
#' `h` (the handling time).
#'
#' The Type III model is defined as:
#' \deqn{ \text{rate} = \frac{b \cdot \text{prey}^2}{1 + b \cdot h \cdot \text{prey}^2} }
#' where \code{b} is the maximum rate of consumption and \code{h} is the handling time.
#'
#' @return A function that calculates the propensity of prey consumption based on the current
#' number of prey and a set of parameters. The returned function takes two arguments:
#' \describe{
#'   \item{prey}{Numeric representing the current number of prey.}
#'   \item{parameters}{A list containing the parameters for the model. Must include:}
#'     \describe{
#'       \item{b}{A numeric value representing the maximum rate of consumption.}
#'       \item{h}{A numeric value representing the handling time.}
#'     }
#' }
#'
#' @examples
#' model <- model_typeIII()
#' parameters <- list(b = 0.3, h = 0.7)
#' model(10, parameters)  # Computes the rate of prey consumption for 10 prey
#'
#' @export
model_typeIII <- function() {
  function(prey, parameters) {
    b <- parameters$b
    h <- parameters$h
    b * prey^2 / (1 + b * h * prey^2)
  }
}

#' Generate a Generalized Holling Model Function for Prey-Predator Simulation
#'
#' This function returns a model function based on the Generalized Holling model, which describes
#' the rate of prey consumption as a function of the number of prey and three parameters:
#' `b` (the maximum rate of consumption), `h` (the handling time), and `q` (the type of functional response).
#'
#' The Generalized Holling model is defined as:
#' \deqn{ \text{rate} = \frac{b \cdot \text{prey}^{(1 + q)}}{1 + b \cdot h \cdot \text{prey}^{(1 + q)} }}
#' where \code{b} is the maximum rate of consumption, \code{h} is the handling time, and \code{q}
#' modifies the type of functional response.
#'
#' @return A function that calculates the propensity of prey consumption based on the current
#' number of prey and a set of parameters. The returned function takes two arguments:
#' \describe{
#'   \item{prey}{Numeric representing the current number of prey.}
#'   \item{parameters}{A list containing the parameters for the model. Must include:}
#'     \describe{
#'       \item{b}{A numeric value representing the maximum rate of consumption.}
#'       \item{h}{A numeric value representing the handling time.}
#'       \item{q}{A numeric value representing the type of functional response.}
#'     }
#' }
#'
#' @examples
#' model <- model_generalised_holling()
#' parameters <- list(b = 0.5, h = 0.2, q = 1)
#' model(10, parameters)  # Computes the rate of prey consumption for 10 prey
#'
#' @export
model_generalised_holling <- function() {
  function(prey, parameters) {
    b <- parameters$b
    h <- parameters$h
    q <- parameters$q
    b * prey^(1 + q) / (1 + b * h * prey^(1 + q))
  }
}
