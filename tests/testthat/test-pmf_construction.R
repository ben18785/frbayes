test_that("create_approx_counts handles basic input correctly", {
  # Create a simple test dataframe
  simulation_result <- tibble(n_prey_remaining = c(1, 2, 2, 3, 4))

  # Test with n_prey_initial = 4
  result <- create_approx_counts(simulation_result, 4)

  expected <- tibble(
    n_prey_remaining = 0:4,
    n = c(0, 1, 2, 1, 1)  # Expected counts for each value
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
    n = c(0, 0, 0, 0)  # All counts should be zero
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
    n = c(1, 3, 2, 4, 0, 0)  # Counts for each value, including those not in the original data
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
