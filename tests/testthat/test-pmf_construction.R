test_that("create_approx_counts handles basic input correctly", {
  # Create a simple test dataframe
  simulation_result <- tibble(n_prey_remaining = c(1, 2, 2, 3, 4))

  # Test with n_prey_initial = 4
  result <- create_approx_counts(simulation_result, 4)

  expected <- tibble(
    n_prey_remaining = 0:4,
    n = c(0, 1, 2, 1, 1) # Expected counts for each value
  )

  expect_equal(result, expected)
})

test_that("create_approx_counts handles empty input", {
  # Test with an empty dataframe
  simulation_result <- tibble(n_prey_remaining = integer(0))

  # Test with n_prey_initial = 3
  result <- create_approx_counts(simulation_result, 3)

  expected <- tibble(
    n_prey_remaining = 0:3,
    n = c(0, 0, 0, 0) # All counts should be zero
  )

  expect_equal(result, expected)
})

test_that("create_approx_counts handles edge cases", {
  # Create a dataframe with some missing values
  simulation_result <- tibble(n_prey_remaining = c(0, 1, 1, 1, 2, 2, 3, 3, 3, 3))

  # Test with n_prey_initial = 5
  result <- create_approx_counts(simulation_result, 5)

  expected <- tibble(
    n_prey_remaining = 0:5,
    n = c(1, 3, 2, 4, 0, 0) # Counts for each value, including those not in the original data
  )

  expect_equal(result, expected)
})

test_that("create_approx_counts handles larger ranges correctly", {
  # Create a dataframe with values from 0 to 10
  simulation_result <- tibble(n_prey_remaining = 0:10)

  # Test with n_prey_initial = 10
  result <- create_approx_counts(simulation_result, 10)

  expected <- tibble(
    n_prey_remaining = 0:10,
    n = rep(1, 11)
  )

  expect_equal(result, expected)
})


test_that("construct_empirical_pmf_log_df calculates PMF and log probabilities correctly", {
  # Create a simple test dataframe
  simulation_result <- tibble(n_prey_remaining = c(1, 2, 2, 3, 4))

  # Test with n_prey_initial = 4 and alpha = 1
  result <- construct_empirical_pmf_log_df(simulation_result, 4, 1)

  # Expected PMF calculation
  expected_counts <- create_approx_counts(simulation_result, 4)
  p <- 1 / (4 + 1)
  H <- rep(p, 5)
  expected_pmf <- (1 * H + expected_counts$n) / sum(1 * H + expected_counts$n)
  expected_log_prob <- log(expected_pmf)

  expected <- tibble(
    n_prey_remaining = 0:4,
    pmf = expected_pmf,
    log_prob = expected_log_prob
  )

  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("construct_empirical_pmf_log_df validates probabilities sum to 1", {
  # Create a simple test dataframe
  simulation_result <- tibble(n_prey_remaining = c(1, 2, 3, 4))

  # Test with n_prey_initial = 4 and alpha = 1
  result <- construct_empirical_pmf_log_df(simulation_result, 4, 2)

  # Check if the sum of PMF is 1
  expect_equal(sum(result$pmf), 1, tolerance = 1e-6)
})

test_that("create_pmf_log returns correct log probability", {
  # Define a sample simulation result
  simulation_result <- tibble(n_prey_remaining = c(1, 2, 2, 3, 4))

  # Create PMF log function
  pmf_log <- create_pmf_log(simulation_result, 4, 1)

  # Expected PMF calculation
  expected_counts <- create_approx_counts(simulation_result, 4)
  p <- 1 / (4 + 1)
  H <- rep(p, 5)
  expected_pmf <- (1 * H + expected_counts$n) / sum(1 * H + expected_counts$n)
  expected_log_prob <- log(expected_pmf)

  # Check log probability values
  for (i in 0:4) {
    expect_equal(pmf_log(i), expected_log_prob[i + 1], tolerance = 1e-6)
  }
})

test_that("create_pmf_log handles out-of-bounds values", {
  # Define a sample simulation result
  simulation_result <- tibble(n_prey_remaining = c(1, 2, 2, 3, 4))

  # Create PMF log function
  pmf_log <- create_pmf_log(simulation_result, 4, 1)

  # Check out-of-bounds values
  expect_equal(pmf_log(-1), -Inf)
  expect_equal(pmf_log(5), -Inf)
})

test_that("create_pmf_log handles invalid inputs", {
  # Define a sample simulation result
  simulation_result <- tibble(n_prey_remaining = c(1, 2, 2, 3, 4))

  # Create PMF log function
  pmf_log <- create_pmf_log(simulation_result, 4, 1)

  # Check invalid inputs
  expect_equal(pmf_log("a"), -Inf)
  expect_equal(pmf_log(2.5), -Inf)
})

test_that("log_probability_single_prey_initial calculates log probability correctly", {
  # Define sample inputs
  parameters <- list(rate = 0.1) # example parameter for the model
  ns_prey_remaining <- c(0, 1, 2, 3, 4)
  n_prey_initial <- 4
  time_max <- 10
  alpha <- 1

  # Create a mock model function
  model <- model_stochastic_degradation()

  # Calculate log probabilities
  n_replicates <- 10000
  result <- log_probability_single_prey_initial(
    parameters = parameters,
    ns_prey_remaining = ns_prey_remaining,
    n_prey_initial = n_prey_initial,
    model = model,
    time_max = time_max,
    n_replicates = n_replicates,
    alpha = alpha
  )

  # Expected result calculation
  simulation_result <- simulate(n_replicates, n_prey_initial, time_max, model, parameters)
  log_prob_function <- create_pmf_log(simulation_result, n_prey_initial, alpha)

  expected_log_prob <- sum(sapply(ns_prey_remaining, log_prob_function))

  expect_equal(result, expected_log_prob, tolerance = 0.5) # simulation yields different pmfs due to stochasticity
})


# Example data for testing
test_data <- tibble(
  n_prey_initial = c(10, 10, 5, 5),
  n_prey_remaining = c(8, 9, 3, 4)
)

test_that("log_probability handles invalid data formats", {
  expect_error(
    log_probability(parameters = list(), data = list(), model = NULL),
    "`data` must be a dataframe or tibble."
  )

  expect_error(
    log_probability(parameters = list(), data = tibble(a = 1:5), model = NULL),
    "`data` must contain the following columns: n_prey_initial, n_prey_remaining"
  )
})

test_that("log_probability handles invalid time_max values", {
  expect_error(
    log_probability(parameters = list(), data = test_data, model = NULL, time_max = -1),
    "`time_max` must be a positive number."
  )

  expect_error(
    log_probability(parameters = list(), data = test_data, model = NULL, time_max = "one"),
    "`time_max` must be a positive number."
  )
})

test_that("log_probability handles invalid n_replicates values", {
  expect_error(
    log_probability(parameters = list(), data = test_data, model = NULL, n_replicates = 0),
    "`n_replicates` must be a positive integer."
  )

  expect_error(
    log_probability(parameters = list(), data = test_data, model = NULL, n_replicates = 1.5),
    "`n_replicates` must be a positive integer."
  )
})

test_that("log_probability handles invalid alpha values", {
  expect_error(
    log_probability(parameters = list(), data = test_data, model = NULL, alpha = -1),
    "`alpha` must be a positive number."
  )

  expect_error(
    log_probability(parameters = list(), data = test_data, model = NULL, alpha = "one"),
    "`alpha` must be a positive number."
  )
})

test_that("log_probability returns the sum of log_probabilities from log_probability_single_prey_initial", {
  # Define parameters for the test
  parameters <- list(rate = 0.1)
  model <- model_stochastic_degradation()
  time_max <- 1
  n_replicates <- 10000
  alpha <- 1

  # Calculate expected log probability manually
  unique_prey_initial <- sort(unique(test_data$n_prey_initial))

  expected_log_prob <- 0
  for (i in seq_along(unique_prey_initial)) {
    df_single_prey_initial <- test_data %>%
      dplyr::filter(n_prey_initial == unique_prey_initial[i])

    ns_prey_remaining <- df_single_prey_initial$n_prey_remaining
    n_prey_initial <- unique_prey_initial[i]

    log_prob_increment <- log_probability_single_prey_initial(
      parameters,
      ns_prey_remaining, n_prey_initial,
      model, time_max,
      n_replicates, alpha
    )

    expected_log_prob <- expected_log_prob + log_prob_increment
  }

  # Get the result from the log_probability function
  result <- log_probability(
    parameters = parameters,
    data = test_data,
    model = model,
    time_max = time_max,
    n_replicates = n_replicates,
    alpha = alpha
  )

  # Check if the result matches the expected log probability
  expect_equal(result, expected_log_prob, tolerance = 0.5)
})
