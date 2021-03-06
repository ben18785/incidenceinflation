#' Samples a case count arising on a given onset day
#'
#' The distribution being drawn from is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' @param max_cases maximum possible cases thought to arise on a given day
#' @inheritParams conditional_cases_logp
#' @param ndraws number of draws of cases
#' @param maximise rather than sample a case count give the case count with the
#' maximum probability (by default is FALSE)
#'
#' @return a sampled case count arising on a given onset day
sample_true_cases_single_onset <- function(
  observation_df, cases_history, max_cases,
  Rt, day_onset, serial_parameters, reporting_parameters,
  ndraws=1, maximise=FALSE) {
  max_observed_cases <- max(observation_df$cases_reported)
  if(max_observed_cases > max_cases)
    stop("Max possible cases should be (much) greater than max observed cases.")
  possible_cases <- max_observed_cases:max_cases
  logps <- conditional_cases_logp(possible_cases, observation_df, cases_history,
                                  Rt, day_onset, serial_parameters, reporting_parameters)
  probs <- exp(logps - matrixStats::logSumExp(logps))
  if(dplyr::last(probs) > 0.01)
    warning(paste0("Cases too few for onset day: ", day_onset,
                   ". Increase max_cases."))
  if(maximise)
    possible_cases[which.max(probs)]
  else
    sample(possible_cases, ndraws, prob=probs, replace=TRUE)
}

#' Calculates max number of days we are uncertain about reporting
#'
#' @param p_gamma_cutoff a p value (0 <= p <= 1) indicating the threshold above which
#' we deem certainty
#' @inheritParams conditional_cases_logp
#'
#' @return a number of days
max_uncertain_days <- function(p_gamma_cutoff, reporting_parameters) {
  r_mean <- reporting_parameters$mean
  r_sd <- reporting_parameters$sd
  days_from_end <- qgamma_mean_sd(p_gamma_cutoff, r_mean, r_sd)
  days_from_end
}

#' Draws a possible history (or histories) of cases
#'
#' The distribution being drawn from at each time t is given by:
#' \deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
#'  p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
#'
#' @param observation_onset_df a tibble with three columns: time_onset, time_reported, cases_reported
#' @inheritParams sample_true_cases_single_onset
#' @inheritParams true_cases
#' @inheritParams max_uncertain_days
#'
#' @return a tibble with an extra cases_estimated column
#' @export
#' @importFrom rlang .data
sample_cases_history <- function(
  observation_onset_df, max_cases,
  Rt_function, serial_parameters, reporting_parameters,
  p_gamma_cutoff=0.99,
  maximise=FALSE) {

  uncertain_period <- max_uncertain_days(p_gamma_cutoff, reporting_parameters)
  start_uncertain_period <- max(observation_onset_df$time_onset) - uncertain_period
  observation_history_df <- observation_onset_df %>%
    dplyr::group_by(.data$time_onset) %>%
    dplyr::mutate(cases_estimated=ifelse(.data$time_onset < start_uncertain_period,
                                    max(.data$cases_reported), NA)) %>%
    dplyr::ungroup()
  onset_times <- unique(observation_history_df$time_onset)
  onset_times_uncertain_period <- onset_times[onset_times >= start_uncertain_period]

  for(i in seq_along(onset_times_uncertain_period)) {
    onset_time <- onset_times_uncertain_period[i]
    snapshots_at_onset_time_df <- observation_history_df %>%
      dplyr::filter(.data$time_onset==onset_time) %>%
      dplyr::select(.data$time_reported, .data$cases_reported)
    pre_observation_df <- observation_history_df %>%
      dplyr::filter(.data$time_onset < onset_time) %>%
      dplyr::select(.data$time_onset, .data$cases_estimated) %>%
      unique() %>%
      dplyr::arrange(dplyr::desc(.data$time_onset))
    cases_history <- pre_observation_df$cases_estimated
    Rt <- Rt_function(onset_time)
    case <- sample_true_cases_single_onset(
      observation_df=snapshots_at_onset_time_df,
      cases_history=cases_history,
      max_cases=max_cases,
      Rt=Rt,
      day_onset=onset_time,
      serial_parameters=serial_parameters,
      reporting_parameters=reporting_parameters,
      ndraws=1,
      maximise=maximise)
    index_onset_time <- which(observation_history_df$time_onset==onset_time)
    observation_history_df$cases_estimated[index_onset_time] <- case
  }
  observation_history_df
}

#' Draws from the gamma distribution or returns the value which maximises
#' it
#'
#' @param shape the shape parameter of a gamma distribution
#' @param rate the rate parameter of a gamma distribution
#' @param ndraws number of draws if maximise=FALSE
#' @param maximise whether to return the mode of the gamma distribution
#'
#' @return a value or (if ndraws > 1) a vector of values
sample_or_maximise_gamma <- function(shape, rate, ndraws, maximise=FALSE) {
  if(maximise)
    (shape - 1) / rate
  else
    stats::rgamma(ndraws, shape, rate)
}

#' Sample a single Rt value corresponding to a single piecewise-
#' constant element of an Rt vector
#'
#' @param Rt_piece_index the index of the Rt piece being sampled
#' @param cases_history_df a tibble with three columns: time_onset, cases_true
#' and Rt_index
#' @param Rt_prior_parameters a list with elements 'shape' and 'rate' describing
#' the gamma prior for Rt
#' @inheritParams sample_cases_history
#' @param ndraws number of draws of Rt
#' @inheritParams true_cases
#'
#' @return a draw (or draws) for Rt
#' @importFrom rlang .data
sample_Rt_single_piece <- function(
    Rt_piece_index, cases_history_df,
    Rt_prior_parameters, serial_parameters,
    serial_max=40, ndraws=1,
    maximise=FALSE) {
  short_df <- cases_history_df %>%
    dplyr::filter(.data$Rt_index <= Rt_piece_index)
  time_max_post_initial_period <- max(short_df$time_onset) - serial_max

  # sample from prior if no data
  posterior_shape <- Rt_prior_parameters$shape
  posterior_rate <- Rt_prior_parameters$rate
  if(time_max_post_initial_period <= 0) {
    return(sample_or_maximise_gamma(
      posterior_shape, posterior_rate, ndraws, maximise))
  }

  # if some data but not enough for whole period
  # do not use truncated points as observed data
  # (but they will be used as covariates)
  short_df <- short_df %>%
    dplyr::mutate(time_after_start = .data$time_onset - serial_max) %>%
    dplyr::mutate(is_observed_data=dplyr::if_else(
      (.data$time_after_start > 0) & (.data$Rt_index == Rt_piece_index), 1, 0))
  onset_times <- short_df %>%
    dplyr::filter(.data$is_observed_data == 1) %>%
    dplyr::pull(.data$time_onset)

  w <- weights_series(serial_max, serial_parameters)
  for(i in seq_along(onset_times)) {
    onset_time <- onset_times[i]
    true_cases <- short_df %>%
      dplyr::filter(.data$time_onset == onset_time) %>%
      dplyr::pull(.data$cases_true)
    posterior_shape <- posterior_shape + true_cases
    cases_history <- short_df %>%
      dplyr::filter(.data$time_onset < onset_time) %>%
      dplyr::arrange(dplyr::desc(.data$time_onset)) %>%
      dplyr::pull(.data$cases_true)
    cases_history <- cases_history[1:serial_max]
    posterior_rate <- posterior_rate + sum(w * cases_history)
  }
  sample_or_maximise_gamma(
    posterior_shape, posterior_rate, ndraws, maximise)
}


#' Sample piecewise-constant Rt values
#'
#' Models the renewal process as from a Poisson:
#' \deqn{cases_true_t ~ Poisson(Rt * \sum_tau=1^t_max w_t cases_true_t-tau))}
#' If an Rt value is given a gamma prior, this results in a posterior
#' distribution:
#' \deqn{Rt ~ gamma(alpha + cases_true_t, beta + \sum_tau=1^t_max w_t cases_true_t-tau))}
#' where alpha and beta are the shape and rate parameters of the gamma
#' prior distribution. Here, we assume that Rt is constant over a set of onset
#' times 'onset_time_set'. This means that the posterior for a single Rt value is given
#' by:
#' \deqn{Rt ~ gamma(alpha + \sum_{t in onset_time_set} cases_true_t,
#'       beta + \sum_{t in onset_time_set}\sum_tau=1^t_max w_t cases_true_t-tau))}
#' This function either returns a draw (or draws if ndraws>1) from
#' this posterior, or it returns the Rt set that maximises it
#' (if maximise=TRUE).
#'
#' @inheritParams sample_Rt_single_piece
#' @return a tibble with three columns: "Rt_piece_index", "draw_index", "Rt"
#' @export
sample_Rt <- function(cases_history_df,
                      Rt_prior_parameters,
                      serial_parameters,
                      serial_max=40,
                      ndraws=1,
                      maximise=FALSE) {
  Rt_piece_indices <- unique(cases_history_df$Rt_index)
  num_Rt_pieces <- length(Rt_piece_indices)
  if(maximise)
    ndraws <- 1
  draw_indices <- seq(1, ndraws, 1)
  m_draws <- matrix(nrow = num_Rt_pieces * ndraws,
                    ncol = 3)
  k <- 1
  for(i in seq_along(Rt_piece_indices)) {
    Rt_piece_index <- Rt_piece_indices[i]
    Rt_vals <- sample_Rt_single_piece(
      Rt_piece_index, cases_history_df,
      Rt_prior_parameters, serial_parameters,
      serial_max, ndraws, maximise=maximise)
    for(j in 1:ndraws) {
      m_draws[k, ] <- c(Rt_piece_index, j, Rt_vals[j])
      k <- k + 1
    }
  }
  colnames(m_draws) <- c("Rt_index", "draw_index", "Rt")
  m_draws <- m_draws %>%
    dplyr::as_tibble()
  m_draws
}

#' Propose new reporting parameters using normal kernel
#' centered at current values
#'
#' @param current_reporting_parameters named list of 'mean' and 'sd' of gamma distribution
#' characterising the reporting delay distribution
#' @param metropolis_parameters named list of 'mean_step', 'sd_step' containing
#' step sizes for Metropolis step
#'
#' @return list of reporting parameters
propose_reporting_parameters <- function(
  current_reporting_parameters,
  metropolis_parameters) {
  mean_now <- current_reporting_parameters$mean
  sd_now <- current_reporting_parameters$sd
  mean_stepsize <- metropolis_parameters$mean_step
  sd_stepsize <- metropolis_parameters$sd_step
  mean_proposed <- stats::rnorm(1, mean_now, mean_stepsize)
  sd_proposed <- stats::rnorm(1, sd_now, sd_stepsize)
  list(mean=mean_proposed, sd=sd_proposed)
}

#' Sample reporting parameters using a single Metropolis step
#'
#' @inheritParams observation_process_all_times_logp
#' @inheritParams propose_reporting_parameters
#'
#' @return list of reporting parameters
metropolis_step <- function(snapshot_with_true_cases_df,
                            current_reporting_parameters,
                            metropolis_parameters) {
  proposed_reporting_parameters <- propose_reporting_parameters(
    current_reporting_parameters,
    metropolis_parameters)
  logp_current <- observation_process_all_times_logp(
    snapshot_with_true_cases_df=snapshot_with_true_cases_df,
    reporting_parameters=current_reporting_parameters
  )
  logp_proposed <- observation_process_all_times_logp(
    snapshot_with_true_cases_df=snapshot_with_true_cases_df,
    reporting_parameters=proposed_reporting_parameters
  )

  log_r <- logp_proposed - logp_current
  log_u <- log(stats::runif(1))
  # nocov start
  if(log_r > log_u)
    proposed_reporting_parameters
  else
    current_reporting_parameters
  # nocov end
}

#' Sample reporting parameters using Metropolis MCMC
#'
#' @inheritParams metropolis_step
#' @param ndraws number of iterates of the Markov chain to simulate
#'
#' @return a tibble with three columns: "draw_index", "mean", "sd"
metropolis_steps <- function(
  snapshot_with_true_cases_df,
  current_reporting_parameters,
  metropolis_parameters,
  ndraws) {

  m_reporting <- matrix(ncol = 3, nrow = ndraws)
  reporting_parameters <- current_reporting_parameters
  for(i in 1:ndraws) {
    reporting_parameters <- metropolis_step(
      snapshot_with_true_cases_df,
      reporting_parameters,
      metropolis_parameters
    )
    m_reporting[i, ] <- c(i,
                          reporting_parameters$mean,
                          reporting_parameters$sd)
  }
  colnames(m_reporting) <- c("draw_index", "mean", "sd")
  m_reporting <- m_reporting %>%
    dplyr::as_tibble()
  m_reporting
}

#' Select reporting parameters by maximising log-probability
#'
#' @inheritParams metropolis_step
#'
#' @return a tibble with three columns: "draw_index", "mean, "sd"
maximise_reporting_logp <- function(
  snapshot_with_true_cases_df,
  current_reporting_parameters) {

  objective_function <- function(theta) {
    -observation_process_all_times_logp(
      snapshot_with_true_cases_df,
      list(mean=theta[1], sd=theta[2]))
  }

  start_point <- c(current_reporting_parameters$mean,
                   current_reporting_parameters$sd)
  theta <- stats::optim(start_point, objective_function)$par
  reporting_parameters <- list(mean=theta[1],
                               sd=theta[2])
  dplyr::tibble(draw_index=1,
         mean=reporting_parameters$mean,
         sd=reporting_parameters$sd)
}

#' Draw reporting parameter values either by sampling or by
#' maximising
#'
#' @inheritParams metropolis_steps
#' @param maximise if true choose reporting parameters by maximising
#' log-probability; else (default) use Metropolis MCMC
#' to draw parameters
#'
#' @return a tibble with three columns: "draw_index", "mean, "sd"
#' @export
sample_reporting <- function(
  snapshot_with_true_cases_df,
  current_reporting_parameters,
  metropolis_parameters,
  maximise=FALSE,
  ndraws=1) {
  if(maximise)
    reporting_parameters <- maximise_reporting_logp(
      snapshot_with_true_cases_df,
      current_reporting_parameters)
  else
    reporting_parameters <- metropolis_steps(
      snapshot_with_true_cases_df,
      current_reporting_parameters,
      metropolis_parameters,
      ndraws=ndraws)

  reporting_parameters
}
