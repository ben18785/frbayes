
create_approx_counts <- function(simulation_result, n_prey_initial) {

  # Count occurrences of each value of n_prey_remaining
    df_counts <- simulation_result %>%
      dplyr::count(n_prey_remaining)

  # Create a tibble with all possible values of n_prey_remaining
  df_all <- dplyr::tibble(n_prey_remaining = 0:n_prey_initial)

  # Merge counts with the tibble of all possible values
  df_all %>%
    dplyr::left_join(df_counts, by = "n_prey_remaining") %>%
    tidyr::replace_na(list(n = 0)) %>%
    dplyr::arrange(n_prey_remaining)
}


construct_empirical_pmf_log_df <- function(simulation_result, n_prey_initial, alpha) {

  # prior measure is uniform on set of possible n_prey_remaining
  p <- 1 / (n_prey_initial + 1)
  H <- rep(p, (n_prey_initial + 1))

  # use Dirichlet Process posterior mean as posterior measure
  df_pmf <- create_approx_counts(simulation_result, n_prey_initial) %>%
    dplyr::mutate(H=H) %>%
    dplyr::mutate(pmf=alpha*H+n) %>%
    dplyr::mutate(pmf=pmf/sum(pmf)) %>%
    dplyr::select(-c("n", "H")) %>%
    dplyr::mutate(log_prob=log(pmf))

  df_pmf
}

create_pmf_log <- function(simulation_result, n_prey_initial, alpha) {

  df_pmf <- construct_empirical_pmf_df(simulation_result, n_prey_initial, alpha)

  pmf_log <- function(n_prey_remaining) {

    if(n_prey_remaining < 0 || n_prey_remaining > n_prey_initial)
      return(-Inf)

    if (!is.numeric(n_prey_remaining) || round(n_prey_remaining) != n_prey_remaining) {
      return(-Inf)
    }

    df_pmf$log_prob[n_prey_remaining + 1]
  }

  pmf_log
}

log_probability_single_prey_initial <- function(
    parameters,
    ns_prey_remaining, n_prey_initial,
    model, time_max,
    n_replicates=1000, alpha=1) {

  simulation_result <- simulate(n_replicates, n_prey_initial, time_max, model, parameters)
  log_prob_function <- create_pmf_log(simulation_result, n_prey_initial, alpha)

  log_prob <- 0
  for(i in seq_along(ns_prey_remaining)) {
    log_prob = log_prob + log_prob_function(ns_prey_remaining[i])
  }

  log_prob
}


log_probability <- function(
    parameters
    data,
    model,
    time_max=1, n_replicates=1000, alpha=1) {

  log_prob <- 0
  unique_prey_initial <- sort(unique(data$n_prey_initial))
  for(i in seq_along(unique_prey_initial)) {

    df_single_prey_initial <- data %>%
      dplyr::filter(n_prey_initial==unique_prey_initial[i])

    log_prob_increment <- log_probability_single_prey_initial(
      parameters,
      ns_prey_remaining, n_prey_initial,
      model, time_max,
      n_replicates, alpha
    )

    log_prob <- log_prob + log_prob_increment

  }

  log_prob
}
