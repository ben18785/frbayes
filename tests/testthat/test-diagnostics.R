test_that("create_bootstrapped_samples handles inputs and outputs correctly", {
  # Mock input data
  data <- data.frame(n_prey_initial = c(10, 20, 30))
  n_bootstraps <- 5
  time_max <- 10
  model <- function(prey, params) {
    1
  }
  mle_parameters <- list(rate = 0.1)

  # Run function
  result <- create_bootstrapped_samples(
    n_bootstraps = n_bootstraps,
    data = data,
    time_max = time_max,
    model = model,
    mle_parameters = mle_parameters
  )

  # Check that the result is a data frame
  expect_s3_class(result, "data.frame")

  # Check that the result contains the expected columns
  expect_true(all(c("n_prey_initial", "n_prey_eaten", "n_prey_remaining", "bootstrap_id") %in% names(result)))

  # Check that the number of rows is correct
  expected_rows <- nrow(data) * n_bootstraps
  expect_equal(nrow(result), expected_rows)

  # Check that bootstrap_id is assigned correctly
  expect_true(all(result$bootstrap_id %in% 1:n_bootstraps))
  expect_equal(length(unique(result$bootstrap_id)), n_bootstraps)
})


test_that("create_bootstrapped_ecdf_real_simulated handles inputs and outputs correctly", {
  # Test with valid inputs
  n_bootstraps <- 5
  data <- data.frame(n_prey_initial = c(10, 20, 30), n_prey_eaten = c(1, 2, 3))
  time_max <- 10
  model <- function(prey, params) {
    1
  }
  mle_parameters <- list(rate = 0.1)

  result <- create_bootstrapped_ecdf_real_simulated(
    n_bootstraps = n_bootstraps,
    data = data,
    time_max = time_max,
    model = model,
    mle_parameters = mle_parameters
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("n_prey_eaten", "ecdf_sim", "ecdf_real", "n_prey_initial", "bootstrap_id") %in% names(result)))

  # Test input validation
  expect_error(create_bootstrapped_ecdf_real_simulated(
    n_bootstraps = -1, data = data, time_max = time_max, model = model, mle_parameters = mle_parameters
  ), "Parameter 'n_bootstraps' must be a positive integer.")

  incomplete_data <- data.frame(n_prey_initial = c(10, 20, 30))
  expect_error(
    create_bootstrapped_ecdf_real_simulated(n_bootstraps = 10, data = incomplete_data, time_max = 10, model = dummy_model, mle_parameters = list(rate = 0.1)),
    "Data frame must contain columns 'n_prey_initial' and 'n_prey_eaten'."
  )

  expect_error(create_bootstrapped_ecdf_real_simulated(
    n_bootstraps = 5, data = "not_a_dataframe", time_max = time_max, model = model, mle_parameters = mle_parameters
  ), "Parameter 'data' must be a data frame.")

  expect_error(create_bootstrapped_ecdf_real_simulated(
    n_bootstraps = 5, data = data, time_max = -1, model = model, mle_parameters = mle_parameters
  ), "Parameter 'time_max' must be a positive numeric value.")

  expect_error(create_bootstrapped_ecdf_real_simulated(
    n_bootstraps = 5, data = data, time_max = time_max, model = "not_a_function", mle_parameters = mle_parameters
  ), "Parameter 'model' must be a function.")

  expect_error(create_bootstrapped_ecdf_real_simulated(
    n_bootstraps = 5, data = data, time_max = time_max, model = model, mle_parameters = NULL
  ), "Parameter 'mle_parameters' must be a non-null list of named parameters.")
})

test_that("create_bootstrapped_ecdf_real_simulated produces correct results", {
  n_bootstraps <- 3
  data <- data.frame(n_prey_initial = c(5, 10), n_prey_eaten = c(2, 4))
  time_max <- 10
  model <- function(prey, params) {
    1
  }
  mle_parameters <- list(rate = 0.1)

  result <- create_bootstrapped_ecdf_real_simulated(
    n_bootstraps = n_bootstraps,
    data = data,
    time_max = time_max,
    model = model,
    mle_parameters = mle_parameters
  )

  # Check that the output is as expected
  expect_equal(nrow(result), n_bootstraps * (6 + 11))
  expect_true("n_prey_eaten" %in% names(result))
  expect_true("ecdf_sim" %in% names(result))
  expect_true("ecdf_real" %in% names(result))
  expect_true("n_prey_initial" %in% names(result))
  expect_true("bootstrap_id" %in% names(result))
})
