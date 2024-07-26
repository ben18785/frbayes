# Set default parameters for testing
params <- list(rate = 100)

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

test_that("model_constant_rate returns a function that calculates correct propensity", {
  # Create the model function
  model <- model_constant_rate()

  # Define parameters
  parameters <- list(rate = 0.1)

  # Test the model function with different numbers of prey
  expect_equal(model(10, parameters), 1.0)  # 0.1 * 10 = 1.0
  expect_equal(model(5, parameters), 0.5)   # 0.1 * 5 = 0.5
  expect_equal(model(0, parameters), 0.0)   # 0.1 * 0 = 0.0

  # Test with different rate
  parameters <- list(rate = 2)
  expect_equal(model(10, parameters), 20.0) # 2 * 10 = 20.0
  expect_equal(model(5, parameters), 10.0)  # 2 * 5 = 10.0

  # Test with different prey values
  parameters <- list(rate = 0.5)
  expect_equal(model(20, parameters), 10.0) # 0.5 * 20 = 10.0
  expect_equal(model(1, parameters), 0.5)   # 0.5 * 1 = 0.5
})

test_that("model_rogersII returns a function that calculates correct propensity", {
  # Create the model function
  model <- model_rogersII()

  # Define parameters
  parameters <- list(a = 0.2, h = 0.5)

  # Test the model function with different numbers of prey
  expect_equal(model(10, parameters), 0.2 * 10 / (1 + 0.2 * 0.5 * 10))  # Expected value calculation
  expect_equal(model(5, parameters), 0.2 * 5 / (1 + 0.2 * 0.5 * 5))    # Expected value calculation
  expect_equal(model(0, parameters), 0.2 * 0 / (1 + 0.2 * 0.5 * 0))    # Expected value calculation

  # Test with different parameters
  parameters <- list(a = 1, h = 1)
  expect_equal(model(10, parameters), 1 * 10 / (1 + 1 * 1 * 10))  # Expected value calculation
  expect_equal(model(5, parameters), 1 * 5 / (1 + 1 * 1 * 5))    # Expected value calculation

  # Test with different prey values
  parameters <- list(a = 0.5, h = 0.2)
  expect_equal(model(20, parameters), 0.5 * 20 / (1 + 0.5 * 0.2 * 20))  # Expected value calculation
  expect_equal(model(1, parameters), 0.5 * 1 / (1 + 0.5 * 0.2 * 1))    # Expected value calculation
})

test_that("model_typeIII returns a function that calculates correct propensity", {
  # Create the model function
  model <- model_typeIII()

  # Define parameters
  parameters <- list(b = 0.3, h = 0.7)

  # Test the model function with different numbers of prey
  expect_equal(model(10, parameters), 0.3 * 10^2 / (1 + 0.3 * 0.7 * 10^2))  # Expected value calculation
  expect_equal(model(5, parameters), 0.3 * 5^2 / (1 + 0.3 * 0.7 * 5^2))    # Expected value calculation
  expect_equal(model(0, parameters), 0.3 * 0^2 / (1 + 0.3 * 0.7 * 0^2))    # Expected value calculation

  # Test with different parameters
  parameters <- list(b = 1, h = 1)
  expect_equal(model(10, parameters), 1 * 10^2 / (1 + 1 * 1 * 10^2))  # Expected value calculation
  expect_equal(model(5, parameters), 1 * 5^2 / (1 + 1 * 1 * 5^2))    # Expected value calculation

  # Test with different prey values
  parameters <- list(b = 0.5, h = 0.2)
  expect_equal(model(20, parameters), 0.5 * 20^2 / (1 + 0.5 * 0.2 * 20^2))  # Expected value calculation
  expect_equal(model(1, parameters), 0.5 * 1^2 / (1 + 0.5 * 0.2 * 1^2))    # Expected value calculation
})
