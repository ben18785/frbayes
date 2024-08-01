# Set default parameters for testing
params <- list(rate = 100)


# this is the previous version of simulate which is easy to read but slow
simulate_reliable <- function(n_replicates, n_prey_initial, time_max, model, parameters) {

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

test_that("simulate_single returns correct number of prey remaining when time_max is not reached", {
  n_prey_initial <- 10
  time_max <- 100
  set.seed(123)  # Set seed for reproducibility
  result <- simulate_single(n_prey_initial, time_max, model_constant_rate(), params)
  expect_equal(result, 0)
})

test_that("simulate_single handles non-positive propensity by stopping simulation", {

  model_non_positive <- function(n_prey, params) {
    return(0)  # Always returns zero propensity
  }

  n_prey_initial <- 10
  time_max <- 5
  expect_warning(result <- simulate_single(n_prey_initial, time_max, model_non_positive, params))
  expect_equal(result, 10)
})


test_that("simulate_single returns initial number of prey when time_max is zero", {
  n_prey_initial <- 10
  time_max <- 0
  result <- simulate_single(n_prey_initial, time_max, model_constant_rate(), params)
  expect_equal(result, 10)
})

test_that("simulate_single handles edge case of zero initial prey", {
  n_prey_initial <- 0
  time_max <- 10
  result <- simulate_single(n_prey_initial, time_max, model_constant_rate(), params)
  expect_equal(result, 0)
})

test_that("simulate_single handles edge case of propensity being non-positive", {
  model_zero_propensity <- function(n_prey, params) {
    return(-1)  # Always returns negative propensity
  }

  n_prey_initial <- 10
  time_max <- 10
  expect_warning(simulate_single(n_prey_initial, time_max, model_zero_propensity, params))
})

test_that("simulate returns a tibble with the correct structure", {
  n_replicates <- 5
  n_prey_initial <- 10
  time_max <- 5
  set.seed(123)  # Set seed for reproducibility

  result <- simulate(n_replicates, n_prey_initial, time_max, model_constant_rate(), params)

  expect_s3_class(result, "tbl_df")
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("replicate_id", "n_prey_initial", "n_prey_eaten", "n_prey_remaining"))
  expect_equal(nrow(result), n_replicates)
})

test_that("simulate validates inputs correctly", {
  expect_error(simulate("a", 10, 5, model_constant_rate(), params), "n_replicates must be a positive integer")
  expect_error(simulate(5, "b", 5, model_constant_rate(), params), "n_prey_initial must be a positive integer")
  expect_error(simulate(5, 10, "c", model_constant_rate(), params), "time_max must be a positive number")
  expect_error(simulate(5, 10, 5, "not_a_function", params), "model must be a function")
  expect_error(simulate(-5, 10, 5, model_constant_rate(), params), "n_replicates must be a positive integer")
  expect_error(simulate(5, -10, 5, model_constant_rate(), params), "n_prey_initial must be a positive integer")
  expect_error(simulate(5, 10, -5, model_constant_rate(), params), "time_max must be a positive number")
})

test_that("simulate returns expected results for given parameters", {
  n_replicates <- 3
  n_prey_initial <- 10
  time_max <- 10
  set.seed(123)  # Set seed for reproducibility

  result <- simulate(n_replicates, n_prey_initial, time_max, model_constant_rate(), params)

  # Since the simulation uses random numbers, we check for expected tibble structure and values within expected ranges
  expect_equal(result$replicate_id, seq_len(n_replicates))
  expect_equal(result$n_prey_initial, rep(n_prey_initial, n_replicates))
  expect_true(all(result$n_prey_eaten >= 0 & result$n_prey_eaten <= n_prey_initial))
  expect_true(all(result$n_prey_remaining >= 0 & result$n_prey_remaining <= n_prey_initial))
})

test_that("simulate approximates correct analytical pmf with constant rates", {

  params <- list(rate = 0.1)
  n_replicates <- 100000
  n_prey_initial <- 10
  time_max <- 5
  set.seed(123)
  result <- simulate(n_replicates, n_prey_initial, time_max, model_constant_rate(), params)

  df_approx <- result %>%
    dplyr::group_by(n_prey_remaining) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(approx=n/sum(n))

  # from eq. 2.7 here: https://arxiv.org/abs/0704.1908
  analytical_pmf <- function(n_prey_remaining, n_prey_initial, t, rate) {
    exp(-rate * n_prey_remaining * t) * (
      choose(n_prey_initial, n_prey_remaining) *
        (1 - exp(-rate * t))^(n_prey_initial-n_prey_remaining))
  }

  pmf <- vector(length = (n_prey_initial + 1))
  for(i in 0:n_prey_initial)
    pmf[i + 1] <- analytical_pmf(i, n_prey_initial, time_max, params$rate)

  df_true <- dplyr::tibble(
    n_prey_remaining=seq(0, n_prey_initial, 1),
    true=pmf
  )

  df_both <- df_approx %>%
    dplyr::left_join(df_true, by="n_prey_remaining") %>%
    dplyr::mutate(abs_diff=abs(true-approx))

  expect_true(all(df_both$abs_diff < 0.1))

})

test_that("simulate behaves similarly to slower but simpler simulation function", {

  params <- list(b = 1.5, h = 0.1, q=1)
  n_replicates <- 100000
  n_prey_initial <- 10
  time_max <- 1
  set.seed(123)
  model <- model_generalised_holling()
  result <- simulate(n_replicates, n_prey_initial, time_max, model, params)

  df_fast <- result %>%
    dplyr::group_by(n_prey_remaining) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fast=n/sum(n))

  result_slow <- simulate_reliable(n_replicates, n_prey_initial, time_max, model, params)
  df_slow <- result_slow %>%
    dplyr::group_by(n_prey_remaining) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(slow=n/sum(n))

  df_both <- df_fast %>%
    dplyr::left_join(df_slow, by="n_prey_remaining") %>%
    dplyr::mutate(abs_diff=abs(fast-slow))

  expect_true(all(df_both$abs_diff < 0.1))

})

test_that("simulate_study throws errors for invalid inputs", {

  # Create a sample valid data frame
  valid_data <- tibble::tibble(
    n_prey_initial = c(10, 20),
    n_replicates = c(5, 10)
  )

  mock_model <- model_constant_rate()
  parameters <- list(rate=0.1)

  # Test invalid `data` input
  expect_error(simulate_study(NULL, 1, mock_model, list()), "`data` must be a dataframe containing columns 'n_prey_initial' and 'n_replicates'.")
  expect_error(simulate_study(tibble::tibble(a = 1, b = 2), 1, mock_model, list()), "`data` must be a dataframe containing columns 'n_prey_initial' and 'n_replicates'.")

  # Test invalid `model` input
  expect_error(simulate_study(valid_data, 1, NULL, list()), "model must be a function")

  # Test invalid `time_max` input
  expect_error(simulate_study(valid_data, -1, mock_model, list()), "time_max must be a positive number")
  expect_error(simulate_study(valid_data, "not numeric", mock_model, list()), "time_max must be a positive number")

  # Test invalid `n_replicates` input
  expect_error(simulate_study(valid_data, 1, mock_model, "not a list"), "`parameters` must be a list")

  # Test invalid dataframe structure
  invalid_data <- tibble::tibble(
    n_prey_initial = c(10, 20, 20),
    n_replicates = c(5, 10, 10)
  )
  expect_error(simulate_study(invalid_data, 1, mock_model, list()), "Dataframe should have one row per n_prey_initial.")
})
