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
