% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmf_construction.R
\name{log_probability}
\alias{log_probability}
\title{Compute the Log Probability of Observed Data Given a Model}
\usage{
log_probability(
  parameters,
  data,
  model,
  time_max = 1,
  n_replicates = 1000,
  alpha = 1
)
}
\arguments{
\item{parameters}{A named list or vector of parameters for the model.}

\item{data}{A data frame or tibble containing the observed data. It must
include at least the following columns:
\itemize{
\item \code{n_prey_initial}: An integer indicating the initial number
of prey.
\item \code{n_prey_remaining}: An integer indicating the number of prey
remaining at the end of the observation period.
}}

\item{model}{A function that calculates the propensity for prey consumption
based on the current state and parameters.}

\item{time_max}{A positive numeric value specifying the maximum time for
the simulation. Defaults to 1.}

\item{n_replicates}{A positive integer specifying the number of replicates
for the simulation. Defaults to 1000.}

\item{alpha}{A positive numeric value representing the Dirichlet Process
hyperparameter. Defaults to 1.}
}
\value{
A numeric value representing the total log probability of the
observed data given the model and parameters.
}
\description{
This function computes the log probability of observing the given data
under a specified model and parameters. The computation is based on
simulations of prey dynamics, using the Dirichlet Process prior and
posterior methods to estimate the probability mass function.
}
\details{
The function first validates the input arguments, including the
data format and value constraints. It then computes the log probability
for each unique initial prey count in the dataset by simulating the
dynamics and using the Dirichlet Process posterior. The final log
probability is the sum of log probabilities across all unique initial
prey counts.
}
\examples{
# Example usage:
parameters <- list(rate = 0.1)
data <- tibble::tibble(
  n_prey_initial = c(10, 20, 30),
  n_prey_remaining = c(5, 10, 15)
)
model <- function(prey, parameters) {
  parameters$rate * prey
}
log_prob <- log_probability(parameters, data, model, time_max = 1, n_replicates = 1000, alpha = 1)
print(log_prob)

}
