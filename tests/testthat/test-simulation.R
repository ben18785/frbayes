# Set default parameters for testing
params <- list(rate = 100)

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

mock_model <- function(n_prey, parameters) {
  # Simple model: returns a constant propensity
  return(0.1)
}

# Define test cases
test_that("simulate_trajectory handles inputs and outputs correctly", {

  # Test valid inputs
  n_prey_initial <- 5
  time_max <- 10
  parameters <- list(rate = 0.1)

  # Run the simulation
  result <- simulate_trajectory(
    n_prey_initial = n_prey_initial,
    time_max = time_max,
    model = model_constant_rate(),
    parameters = parameters
  )

  # Check that the result is a data frame
  expect_s3_class(result, "data.frame")

  # Check that the result contains the expected columns
  expect_true(all(c("time", "n_prey_remaining", "n_prey_eaten") %in% names(result)))

  # Check that 'n_prey_remaining' decreases over time
  expect_true(all(diff(result$n_prey_remaining) <= 0))
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

  # test simulate runs with n_prey_initial=1
  n_prey_initial <- 1
  result <- simulate(n_replicates, n_prey_initial, time_max, model_constant_rate(), params)
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

test_that("simulate_study produces reasonable data", {

  experimental_setup <- data.frame(
    n_prey_initial = c(5, 10, 20, 30, 40),
    n_replicates = 100
  )

  # generate synthetic data
  true_parameters <- list(a=2, h=0.1)
  df <- simulate_study(
    data=experimental_setup,
    time_max = 1,
    model = model_rogersII(),
    parameters = true_parameters
  )

  expect_equal(nrow(df), length(experimental_setup$n_prey_initial) * experimental_setup$n_replicates[1])
  cnames_exp <- c("replicate_id", "n_prey_initial", "n_prey_eaten", "n_prey_remaining")
  expect_true(all(cnames_exp %in% colnames(df)))
})
