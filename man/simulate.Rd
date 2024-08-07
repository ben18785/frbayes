% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulate}
\alias{simulate}
\title{Simulate Multiple Runs of a Prey-Predator Model}
\usage{
simulate(n_replicates, n_prey_initial, time_max, model, parameters)
}
\arguments{
\item{n_replicates}{An integer representing the number of replicate simulations to perform.}

\item{n_prey_initial}{An integer representing the initial number of prey.}

\item{time_max}{A numeric value indicating the maximum time for the simulation.}

\item{model}{A function that calculates the propensity given the current number of prey and parameters.}

\item{parameters}{A list of named parameters required by the model function.}
}
\value{
A tibble with columns \code{n_prey_initial}, \code{n_prey_eaten}, and \code{n_prey_remaining},
summarizing the results of the simulations. Each row represents a single simulation replicate.
}
\description{
This function simulates multiple runs of a prey-predator model, each with the same
initial conditions and parameters, and returns a summary of the results.
}
\examples{
# Define a simple linear model for the propensity
model <- function(n_prey, params) {
  rate <- params$rate
  return(rate * n_prey)
}

# Set the number of replicates, initial number of prey, maximum time, and model parameters
n_replicates <- 10
n_prey_initial <- 100
time_max <- 10
parameters <- list(rate = 0.1)

# Run the simulations
simulate(n_replicates, n_prey_initial, time_max, model, parameters)

}
