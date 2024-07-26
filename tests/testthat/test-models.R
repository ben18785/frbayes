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
