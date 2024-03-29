test_that("sample_true_cases_single_onset produces reasonable case draws", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=5, sd=3)
  max_cases <- 100
  s_params <- list(mean=10, sd=1)
  Rt <- 2
  t_max <- 30
  cases_history <- rep(4, t_max)
  day_onset <- 0
  w <- weights_series(t_max, s_params)
  mean_cases <- expected_cases(Rt, w, cases_history)
  case <- sample_true_cases_single_onset(observation_df=observation_matrix,
                                 cases_history=cases_history,
                                 max_cases=max_cases,
                                 Rt=Rt,
                                 day_onset=day_onset,
                                 serial_parameters=s_params,
                                 reporting_parameters=reporting_parameters)
  expect_true(case >= min(observation_matrix$cases_reported))
  expect_true(case <= max_cases)

  # when reporting mean is low, expect to have pretty much seen all cases
  reporting_parameters <- list(mean=1, sd=1)
  cases <- sample_true_cases_single_onset(observation_df=observation_matrix,
                            cases_history=cases_history,
                            max_cases=max_cases,
                            Rt=Rt,
                            day_onset=day_onset,
                            serial_parameters=s_params,
                            reporting_parameters=reporting_parameters,
                            ndraws=20)
  expect_true(any(cases <= 5))

  # check that cases jump if serial interval looks further back
  cases_history <- c(rep(4, t_max / 2), rep(1000, t_max / 2))
  s_params <- list(mean=1, sd=1)
  cases <- sample_true_cases_single_onset(observation_df=observation_matrix,
                             cases_history=cases_history,
                             max_cases=max_cases,
                             Rt=Rt,
                             day_onset=day_onset,
                             serial_parameters=s_params,
                             reporting_parameters=reporting_parameters,
                             ndraws=20)
  s_params <- list(mean=25, sd=1)
  cases1 <- sample_true_cases_single_onset(observation_df=observation_matrix,
                             cases_history=cases_history,
                             max_cases=max_cases,
                             Rt=Rt,
                             day_onset=day_onset,
                             serial_parameters=s_params,
                             reporting_parameters=reporting_parameters,
                             ndraws=20)
  expect_true(mean(cases1) > mean(cases))

  max_cases <- 2
  expect_warning(sample_true_cases_single_onset(observation_df=observation_matrix,
                                                cases_history=cases_history,
                                                max_cases=max_cases,
                                                Rt=Rt,
                                                day_onset=day_onset,
                                                serial_parameters=s_params,
                                                reporting_parameters=reporting_parameters,
                                                ndraws=20))
  max_cases <- 0
  expect_error(sample_true_cases_single_onset(
    observation_df=observation_matrix,
    cases_history=cases_history,
    max_cases=max_cases,
    Rt=Rt,
    day_onset=day_onset,
    serial_parameters=s_params,
    reporting_parameters=reporting_parameters,
    ndraws=20))
})

test_that("sample_true_cases_single_onset produces consistent maximum", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=5, sd=3)
  max_cases <- 100
  s_params <- list(mean=10, sd=1)
  Rt <- 2
  t_max <- 30
  cases_history <- rep(4, t_max)
  day_onset <- 0
  max_observed_cases <- max(observation_matrix$cases_reported)
  possible_cases <- max_observed_cases:max_cases
  logps <- conditional_cases_logp(possible_cases, observation_matrix, cases_history,
                                  Rt, day_onset, s_params, reporting_parameters)

  map_val <- possible_cases[which.max(logps)]
  case <- sample_true_cases_single_onset(
    observation_df=observation_matrix,
    cases_history=cases_history,
    max_cases=max_cases,
    Rt=Rt,
    day_onset=day_onset,
    serial_parameters=s_params,
    reporting_parameters=reporting_parameters,
    maximise=T)
  expect_equal(case, map_val)
})

test_that("max_uncertain_days selects reasonable gamma distribution quanties", {
  r_params <- list(mean=7, sd=0.01)
  expect_equal(round(max_uncertain_days(0.5, r_params)), 7)
  expect_true(max_uncertain_days(0.95, r_params) > 7)

  r_params <- list(mean=7, sd=5)
  expect_true(max_uncertain_days(0.95, r_params) > 10)
  expect_equal(max_uncertain_days(0, r_params), 0)
})

test_that("sample_cases_history adds cases_estimated that look reasonable", {
  days_total <- 100
  kappa <- 1000
  r_params <- list(mean=10, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T) %>%
    dplyr::mutate(reporting_piece_index=1)
  max_cases <- 5000
  r_params <- dplyr::tibble(mean=10, sd=3, reporting_piece_index=1)
  df_est <- sample_cases_history(df, max_cases, Rt_function, s_params, r_params)
  expect_true(all.equal(df_est %>% dplyr::select(-cases_estimated),
                        df))
  expect_true(max(df_est$cases_estimated) < max_cases)
  expect_true(min(df_est$cases_estimated) > 0)
  expect_equal(sum(is.na(df_est$cases_estimated)), 0)

  df_group <- df_est %>%
    dplyr::group_by(time_onset) %>%
    dplyr::summarise(
      cases_reported=max(cases_reported),
      cases_estimated=min(cases_estimated),
      cases_true=mean(cases_true)) %>%
    dplyr::mutate(diff=cases_estimated-cases_reported) %>%
    dplyr::mutate(diff_true=cases_estimated-cases_true)
  expect_equal(sum(df_group$diff < 0), 0)
  expect_true(max(abs(df_group$diff_true)) < 300)

  # throws error when reporting_piece_index not in observation_onset_df
  df_tmp <- df %>%
    dplyr::select(-reporting_piece_index)
  expect_error(sample_cases_history(df_tmp, max_cases, Rt_function, s_params, r_params))
})

test_that("sample_cases_history yields a single case history when
maximising", {
  days_total <- 30
  kappa <- 1000
  r_params <- list(mean=10, sd=3)
  v_Rt <- c(rep(1.5, 10), rep(0.4, 10), rep(1.5, 10))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  s_params <- list(mean=5, sd=3)
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T) %>%
    dplyr::mutate(reporting_piece_index=1)
  max_cases <- 5000
  r_params <- dplyr::tibble(
    reporting_piece_index=1,
    mean=r_params$mean,
    sd=r_params$sd
  )
  f_est <- function(i) {
    df_est <- sample_cases_history(df, max_cases, Rt_function,
                         s_params, r_params,
                         maximise = T)
    df_est$cases_estimated
  }
  cases <- purrr::map(seq(1, 4, 1), f_est)
  expect_true(all.equal(cases[[1]], cases[[2]]))
  expect_true(all.equal(cases[[2]], cases[[3]]))
  expect_true(all.equal(cases[[3]], cases[[4]]))
})

test_that("sample_or_maximise_from_gamma returns either draws or maximum of
          gamma", {
  shape <- 5
  rate <- 5
  ndraws <- 1
  val <- sample_or_maximise_gamma(shape, rate, ndraws)
  expect_equal(length(val), ndraws)
  ndraws <- 20
  vals <- sample_or_maximise_gamma(shape, rate, ndraws)
  expect_equal(length(vals), ndraws)

  # maximise
  val <- sample_or_maximise_gamma(shape, rate, ndraws,
                                  maximise=T)
  expect_equal(val, (shape - 1) / rate)
})

# tests for sampling Rt
days_total <- 100
Rt_1 <- 1.5
Rt_2 <- 1.0
Rt_3 <- 1.3
v_Rt <- c(rep(Rt_1, 40), rep(Rt_2, 20), rep(Rt_3, 40))
Rt_function <- stats::approxfun(1:days_total, v_Rt)
s_params <- list(mean=5, sd=3)
r_params <- list(mean=10, sd=3)
kappa <- 1000
df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                         kappa=kappa) %>%
  dplyr::select(time_onset, cases_true)
Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))
Rt_index_lookup <- dplyr::tibble(
  time_onset=seq_along(Rt_indices),
  Rt_index=Rt_indices)
df <- df %>%
  dplyr::left_join(Rt_index_lookup, by = "time_onset") %>%
  dplyr::select(time_onset, cases_true, Rt_index) %>%
  unique()
Rt_prior <- list(shape=1, rate=1)

test_that("sample_Rt_single_piece returns reasonable values: these are
          basically functional tests", {

  ndraws <- 42
  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = ndraws,
    serial_max = 20)
  expect_equal(ndraws, length(Rt_vals))

  f_Rt <- function(piece) {
    sample_Rt_single_piece(
      piece, df,
      Rt_prior, s_params,
      ndraws = 1000,
      serial_max = 20)
  }

  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 40)

  # test for maximisation
  f_Rt <- function(i) sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 20,
    maximise = T)
  expect_equal(length(f_Rt(1)), 1)
  Rt_vals1 <- purrr::map_dbl(seq(1, 3, 1), f_Rt)
  # no variation in maximum so should be no sd
  expect_equal(sd(Rt_vals1), 0)

  # mode of sampling distribution should be close to max
  mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = 1000,
    serial_max = 20)
  expect_true(abs(mode(Rt_vals) - Rt_vals[1]) < 0.2)
})

test_that("sample_Rt_single_piece works ok with NB renewal model", {
  ndraws <- 420
  Rt_vals <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = ndraws,
    serial_max = 20)

  # for large kappa draws should be similar
  kappa <- 10000
  Rt_vals_nb <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = ndraws,
    serial_max = 20,
    kappa=kappa,
    is_negative_binomial=TRUE)
  expect_true(abs(mean(Rt_vals_nb) - mean(Rt_vals)) < 0.5)
  expect_true(abs(sd(Rt_vals_nb) - sd(Rt_vals)) < 0.1)

  # for small kappa the NB sd should exceed Poisson
  kappa <- 0.5
  Rt_vals_nb <- sample_Rt_single_piece(
    2, df,
    Rt_prior, s_params,
    ndraws = ndraws,
    serial_max = 20,
    kappa=kappa,
    is_negative_binomial=TRUE)
  expect_true(sd(Rt_vals_nb) > sd(Rt_vals))
})

test_that("sample_Rt returns sensible values", {
  ndraws <- 300
  Rt_df <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20)

  # check shapes of outputs
  cnames <- colnames(Rt_df)
  expect_true(all.equal(cnames,
                        c("Rt_index", "draw_index", "Rt")))
  npieces <- max(df$Rt_index)
  expect_equal(npieces * ndraws, nrow(Rt_df))

  # check substance of outputs
  f_Rt <- function(piece) {
    val <- sample_Rt_single_piece(
      piece, df,
      Rt_prior, s_params,
      ndraws = 1000,
      serial_max = 20)
    mean(val)
  }
  indices <- seq(1, 5, 1)
  Rt_vals <- purrr::map_dbl(indices, f_Rt)
  single_df <- dplyr::tibble(Rt_index=indices,
                             Rt_single=Rt_vals)
  Rt_df <- Rt_df %>%
    dplyr::group_by(.data$Rt_index) %>%
    dplyr::summarise(Rt=mean(Rt)) %>%
    dplyr::left_join(single_df, by = "Rt_index") %>%
    dplyr::mutate(diff=Rt-Rt_single)
  expect_true(abs(sum(Rt_df$diff)) < 0.2)

  # check maximisation
  f_Rt <- function(piece) {
    val <- sample_Rt_single_piece(
      piece, df,
      Rt_prior, s_params,
      serial_max = 20,
      maximise = T)
    mean(val)
  }
  indices <- seq(1, 5, 1)
  Rt_vals <- purrr::map_dbl(indices, f_Rt)
  single_df <- dplyr::tibble(Rt_index=indices,
                             Rt_single=Rt_vals)
  Rt_df <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20,
                     maximise = T)
  expect_equal(nrow(Rt_df), max(df$Rt_index))
  Rt_df <- Rt_df %>%
    dplyr::left_join(single_df, by = "Rt_index") %>%
    dplyr::mutate(diff=Rt-Rt_single)
  expect_true(abs(sum(Rt_df$diff)) < 0.0001)
})

test_that("sample_Rt works for negative binomial renewal model", {

  ndraws <- 300
  Rt_df <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20)

  Rt_df_nb <- sample_Rt(df,
                     Rt_prior, s_params,
                     ndraws = ndraws,
                     serial_max = 20,
                     kappa=1,
                     is_negative_binomial=TRUE)
  expect_equal(nrow(Rt_df), nrow(Rt_df_nb))
  expect_equal(ncol(Rt_df), ncol(Rt_df_nb))

  # errors thrown if not providing kappa
  expect_error(sample_Rt(df,
                         Rt_prior, s_params,
                         ndraws = ndraws,
                         serial_max = 20,
                         is_negative_binomial=TRUE))

  # errors thrown if providing invalid kappa
  expect_error(sample_Rt(df,
                         Rt_prior, s_params,
                         ndraws = ndraws,
                         serial_max = 20,
                         kappa=0,
                         is_negative_binomial=TRUE))
})

test_that("accept_reject works ok", {

  logp_current <- 0
  logp_proposed <- -10000
  current_parameters <- 0
  proposed_parameters <- 1
  list_param_logp <- accept_reject(logp_current, logp_proposed,
    current_parameters,
    proposed_parameters)
  expect_equal(current_parameters, list_param_logp$parameter)
  expect_equal(logp_current, list_param_logp$logp)
})

test_that("prior_reporting_parameters works ok", {
  current_reporting_parameters <- list(mean=3, sd=2, reporting_piece_index=1)
  prior_params <- list(mean_mu=3, mean_sigma=2,
                       sd_mu=5, sd_sigma=1)
  val <- prior_reporting_parameters(current_reporting_parameters,
                                    prior_params)
  val1 <- dgamma_mean_sd(current_reporting_parameters$mean,
                         prior_params$mean_mu,
                         prior_params$mean_sigma,
                         log=TRUE)
  val2 <- dgamma_mean_sd(current_reporting_parameters$sd,
                         prior_params$sd_mu,
                         prior_params$sd_sigma,
                         log=TRUE)
  expect_equal(val, val1 + val2)
})

test_that("propose_reporting_parameters works ok", {
  r_params <- dplyr::tibble(mean=4, sd=4, reporting_piece_index=1)
  met_params <- list(mean_step=0.1, sd_step=0.2)
  r_params1 <- propose_reporting_parameters(
    r_params, met_params)
  expect_true(abs(r_params$mean - r_params1$mean) < 3)
  expect_true(abs(r_params$sd - r_params1$sd) < 3)

  met_params <- list(mean_step=0.000001, sd_step=0.2)
  r_params1 <- propose_reporting_parameters(
    r_params, met_params)
  expect_true(abs(r_params$mean - r_params1$mean) < 0.1)
})

days_total <- 30
df <- generate_snapshots(days_total, Rt_function,
                         s_params, r_params,
                         kappa=kappa, thinned=T) %>%
  dplyr::mutate(reporting_piece_index=1)

test_that("metropolis_step works as expected", {
  met_params <- list(mean_step=0.01, sd_step=0.01)
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  r_params <- dplyr::tibble(mean=r_params$mean,
                     sd=r_params$sd,
                     reporting_piece_index=1)
  logp_current <- -1000
  list_param_logp <- metropolis_step(
    df, r_params, logp_current,
    prior_params, met_params)
  r_params1 <- list_param_logp$parameter
  expect_true(abs(r_params1$mean - r_params$mean) < 0.2)
  expect_true(abs(r_params1$sd - r_params$sd) < 0.2)
})

test_that("metropolis_steps returns multiple steps", {

  met_params <- list(mean_step=0.01, sd_step=0.01)
  ndraws <- 10
  rep_prior_params <- list(mean_mu=5, sd_mu=3,
                           mean_sigma=5, sd_sigma=3)
  r_params <- dplyr::tibble(mean=r_params$mean,
                     sd=r_params$sd,
                     reporting_piece_index=1)
  logp_current <- -1000
  list_draws_logp <- metropolis_steps(
    snapshot_with_true_cases_df=df,
    current_reporting_parameters=r_params,
    logp_current=logp_current,
    prior_parameters=rep_prior_params,
    metropolis_parameters=met_params,
    ndraws=ndraws)
  expect_equal(names(list_draws_logp)[2], "logp")
  output <- list_draws_logp$reporting_parameters
  expect_equal(nrow(output), ndraws)
  expect_true(all.equal(colnames(output),
                        c("reporting_piece_index", "draw_index", "mean", "sd")))
})

test_that("propose_overdispersion_parameter works ok", {
  current <- 0
  sd <- 2
  new <- propose_overdispersion_parameter(current, sd)
  expect_true(current != new)
})

test_that("metropolis_step_overdispersion works ok", {

  overdispersion_current <- 100
  logp_current <- -100000
  cases_history_rt_df <- df %>%
    dplyr::select("time_onset", "cases_true") %>%
    unique() %>%
    dplyr::mutate(Rt=2)

  list_param_logp <- metropolis_step_overdispersion(
    overdispersion_current,
    logp_current,
    cases_history_rt_df,
    serial_parameters=list(mean=2, sd=3),
    prior_overdispersion_parameter=list(mean=10, sd=100),
    overdispersion_metropolis_sd=1)
  expect_true(list_param_logp$overdispersion != overdispersion_current)
  expect_true(list_param_logp$logp > logp_current)
})

test_that("maximise_reporting_logp maximises prob", {
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  r_params <- dplyr::tibble(reporting_piece_index=1,
                            mean=r_params$mean,
                            sd=r_params$sd)
  df <- df %>%
    dplyr::mutate(reporting_piece_index=1)
  output <- maximise_reporting_logp(df, r_params, prior_params)
  expect_true(abs(output$mean - r_params$mean) < 0.4)
  expect_true(abs(output$sd - r_params$sd) < 0.4)
  expect_equal(nrow(output), 1)
})

test_that("maximise_reporting_logp throws errors with piece info missing", {

  df <- df %>%
    dplyr::select(-"reporting_piece_index")
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)

  # r_params doesn't contain reporting_piece_index
  df_temp <- df %>%
    dplyr::mutate(reporting_piece_index=1)
  expect_error(maximise_reporting_logp(df_temp, r_params, prior_params))

  # df doesn't contain reporting_piece_index
  r_params <- dplyr::tibble(reporting_piece_index=1,
                            mean=r_params$mean,
                            sd=r_params$sd)
  expect_error(maximise_reporting_logp(df, r_params, prior_params))

})

test_that("maximise_reporting_logp works ok with multiple pieces", {

  days_total <- 30
  r_params <- dplyr::tibble(mean=c(rep(10, 15), rep(3, 15)),
                     sd=c(rep(3, 15), rep(1, 15)),
                     time_onset=seq_along(mean)) %>%
    dplyr::mutate(reporting_piece_index=c(rep(1, 15), rep(2, 15)))
  df <- generate_snapshots(days_total, Rt_function,
                           s_params, r_params,
                           kappa=kappa, thinned=T)
  r_params_short <- r_params %>%
    dplyr::select(time_onset, reporting_piece_index) %>%
    unique()

  df <- df %>%
    dplyr::left_join(r_params_short, by="time_onset")
  r_params <- dplyr::tibble(reporting_piece_index=c(1, 2),
                     mean=c(8, 4),
                     sd=c(4, 2))
  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  output <- maximise_reporting_logp(df, r_params, prior_params)

  # try optimisation on each subperiod: should yield same
  r_params_1 <- r_params %>%
    dplyr::filter(reporting_piece_index==1)
  df_1 <- df %>%
    dplyr::filter(reporting_piece_index==1)
  output_1 <- maximise_reporting_logp(df_1, r_params_1, prior_params)
  r_params_2 <- r_params %>%
    dplyr::filter(reporting_piece_index==2)
  df_2 <- df %>%
    dplyr::filter(reporting_piece_index==2)
  output_2 <- maximise_reporting_logp(df_2, r_params_2, prior_params)

  expect_equal(output$mean[1], output_1$mean[1])
  expect_equal(output$mean[2], output_2$mean[1])
  expect_equal(output$sd[1], output_1$sd[1])
  expect_equal(output$sd[2], output_2$sd[1])
  expect_equal(output$reporting_piece_index[1], output_1$reporting_piece_index[1])
  expect_equal(output$reporting_piece_index[2], output_2$reporting_piece_index[1])
})

test_that("sample_reporting produces output of correct shape", {

  prior_params <- list(mean_mu=2, mean_sigma=100,
                       sd_mu=2, sd_sigma=100)
  met_params <- list(mean_step=0.01, sd_step=0.01)
  r_params <- dplyr::tibble(mean=r_params$mean,
                     sd=r_params$sd,
                     reporting_piece_index=1)
  logp_current <- -1000
  result <- sample_reporting(df, r_params, logp_current, prior_params, met_params)
  expect_equal(names(result)[2], "logp")
  output <- result$reporting_parameters
  expect_equal(nrow(output), 1)
  expect_equal(max(output$draw_index), 1)

  ndraws <- 2
  result <- sample_reporting(df, r_params, logp_current, prior_params, met_params,
                             ndraws=ndraws)
  output <- result$reporting_parameters
  expect_equal(nrow(output), ndraws)
  expect_equal(max(output$draw_index), ndraws)

  result <- sample_reporting(df, r_params, logp_current, prior_params, met_params,
                             maximise=T)
  output <- result$reporting_parameters
  expect_equal(nrow(output), 1)
  expect_equal(max(output$draw_index), 1)
})


test_that("mcmc produces outputs of correct shape", {

  niter <- 2
  days_total <- 100
  r_params <- list(mean=10, sd=3)
  s_params <- list(mean=5, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  Rt_prior <- list(shape=1, rate=1)
  kappa <- 10
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                           kappa=kappa)

  snapshot_with_Rt_index_df <- df
  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
    dplyr::select(-cases_true)
  initial_Rt <- tidyr::tribble(~Rt_index, ~Rt,
                        1, 1.5,
                        2, 1.5,
                        3, 0.4,
                        4, 1.5,
                        5, 1.5)
  Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))

  Rt_index_lookup <- tidyr::tibble(
    time_onset=seq_along(Rt_indices),
    Rt_index=Rt_indices)
  snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
    dplyr::left_join(Rt_index_lookup, by="time_onset") %>%
    dplyr::mutate(reporting_piece_index=1)

  initial_reporting_parameters <- dplyr::tibble(mean=5, sd=3, reporting_piece_index=1)
  serial_parameters <- list(mean=5, sd=3)
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  # test throws errors
  ## wrongly named cols
  wrong_df <- snapshot_with_Rt_index_df %>%
    dplyr::rename(time_onset_wrong=time_onset)
  expect_error(mcmc(niterations=niter,
                    wrong_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE))

  ## too many cols
  wrong_df <- snapshot_with_Rt_index_df %>%
    dplyr::mutate(unnecessary_col="hi")
  expect_error(mcmc(niterations=niter,
                    wrong_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE))

  # five columns but not reporting_piece_index
  wrong_df <- snapshot_with_Rt_index_df %>%
    dplyr::rename(reporting_piece_index_wrong=reporting_piece_index)
  expect_error(mcmc(niterations=niter,
                    wrong_df,
                    priors,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE))

  # overdispersion parameter <= 0
  expect_error(mcmc(niterations=niter,
                   snapshot_with_Rt_index_df,
                   priors,
                   serial_parameters,
                   initial_cases_true,
                   initial_reporting_parameters,
                   initial_Rt,
                   reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                   serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
                   initial_overdispersion = -1))

  # maximisation with a NB model doesn't work
  priors_tmp <- priors
  priors_tmp$overdispersion <- list(mean=5, sd=10)
  expect_error(mcmc(niterations=niter,
                    snapshot_with_Rt_index_df,
                    priors_tmp,
                    serial_parameters,
                    initial_cases_true,
                    initial_reporting_parameters,
                    initial_Rt,
                    reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
                    serial_max=40, p_gamma_cutoff=0.99, maximise=TRUE,
                    is_negative_binomial=TRUE))


  # test MCMC sampling
  niter <- 5
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE)

  ## check outputs
  ### overall
  expect_equal(length(res), 3)

  ### cases
  cases_df <- res$cases
  expect_true(all.equal(c("time_onset", "cases_true", "iteration", "chain"),
                        colnames(cases_df)))
  expect_equal(min(cases_df$iteration), 1)
  expect_equal(max(cases_df$iteration), niter)
  expect_equal(min(cases_df$time_onset), min(df$time_onset))
  expect_equal(max(cases_df$time_onset), max(df$time_onset))
  expect_equal(min(cases_df$chain), 1)
  expect_equal(max(cases_df$chain), 1)

  ## Rt
  rt_df <- res$Rt
  expect_true(all.equal(c("iteration", "Rt_index", "Rt", "chain"),
                        colnames(rt_df)))
  expect_equal(min(rt_df$iteration), 1)
  expect_equal(max(rt_df$iteration), niter)
  expect_equal(min(rt_df$Rt_index), min(initial_Rt$Rt_index))
  expect_equal(max(rt_df$Rt_index), max(initial_Rt$Rt_index))
  expect_equal(min(rt_df$chain), 1)
  expect_equal(max(rt_df$chain), 1)

  # reporting delays
  reporting_df <- res$reporting
  expect_true(all.equal(c("reporting_piece_index", "mean", "sd", "iteration", "chain"),
                        colnames(reporting_df)))
  expect_equal(min(reporting_df$iteration), 1)
  expect_equal(max(reporting_df$iteration), niter)
  expect_equal(min(reporting_df$chain), 1)
  expect_equal(max(reporting_df$chain), 1)

  # test optimisation
  niter <- 2 # needed since maximisation is iterative
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=TRUE)
  ### overall
  expect_equal(length(res), 3)

  ### cases
  cases_df <- res$cases
  expect_true(all.equal(c("time_onset", "cases_true", "iteration", "chain"),
                        colnames(cases_df)))
  expect_equal(min(cases_df$iteration), 1)
  expect_equal(max(cases_df$iteration), niter)
  expect_equal(min(cases_df$time_onset), min(df$time_onset))
  expect_equal(max(cases_df$time_onset), max(df$time_onset))
  expect_equal(min(cases_df$chain), 1)
  expect_equal(max(cases_df$chain), 1)

  ## Rt
  rt_df <- res$Rt
  expect_true(all.equal(c("iteration", "Rt_index", "Rt", "chain"),
                        colnames(rt_df)))
  expect_equal(min(rt_df$iteration), 1)
  expect_equal(max(rt_df$iteration), niter)
  expect_equal(min(rt_df$Rt_index), min(initial_Rt$Rt_index))
  expect_equal(max(rt_df$Rt_index), max(initial_Rt$Rt_index))
  expect_equal(min(rt_df$chain), 1)
  expect_equal(max(rt_df$chain), 1)

  # reporting delays
  reporting_df <- res$reporting
  expect_true(all.equal(c("reporting_piece_index", "mean", "sd", "iteration", "chain"),
                        colnames(reporting_df)))
  expect_equal(min(reporting_df$iteration), 1)
  expect_equal(max(reporting_df$iteration), niter)
  expect_equal(min(reporting_df$chain), 1)
  expect_equal(max(reporting_df$chain), 1)


  # test MCMC sampling with NB model
  niter <- 5

  ## prior not provided on overdispersion parameter
  expect_error(mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              is_negative_binomial=TRUE))

  # sampling works if given overdispersion prior
  priors$overdispersion <- list(mean=5, sd=2)
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              is_negative_binomial=TRUE)
  expect_equal(length(res), 4)
  expect_true("overdispersion" %in% names(res))
  overdispersion <- res$overdispersion
  expect_equal(nrow(overdispersion), 5)
  expect_true("overdispersion" %in% names(overdispersion))
  expect_true("iteration" %in% names(overdispersion))
})

test_that("multiple chains works", {

  days_total <- 100
  r_params <- list(mean=10, sd=3)
  s_params <- list(mean=5, sd=3)
  v_Rt <- c(rep(1.5, 40), rep(0.4, 20), rep(1.5, 40))
  Rt_function <- stats::approxfun(1:days_total, v_Rt)
  Rt_prior <- list(shape=1, rate=1)
  kappa <- 10
  df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                           kappa=kappa)

  snapshot_with_Rt_index_df <- df
  initial_cases_true <- df %>%
    dplyr::select(time_onset, cases_true) %>%
    unique()
  snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
    dplyr::select(-cases_true)
  initial_Rt <- tidyr::tribble(~Rt_index, ~Rt,
                               1, 1.5,
                               2, 1.5,
                               3, 0.4,
                               4, 1.5,
                               5, 1.5)
  Rt_indices <- unlist(purrr::map(seq(1, 5, 1), ~rep(., 20)))

  Rt_index_lookup <- tidyr::tibble(
    time_onset=seq_along(Rt_indices),
    Rt_index=Rt_indices)
  snapshot_with_Rt_index_df <- snapshot_with_Rt_index_df %>%
    dplyr::left_join(Rt_index_lookup)

  initial_reporting_parameters <- list(mean=5, sd=3)
  serial_parameters <- list(mean=5, sd=3)
  priors <- list(Rt=Rt_prior,
                 reporting=list(mean_mu=5,
                                mean_sigma=10,
                                sd_mu=3,
                                sd_sigma=5),
                 max_cases=5000)

  # multiple chains in serial
  niter <- 2
  nchains <- 2
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              nchains=nchains)
  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)

  nchains <- 3
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              nchains=nchains)
  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)

  # multiple chains in parallel
  library(doParallel)
  cl <- makeCluster(2)
  registerDoParallel(cl)
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              nchains=nchains, is_parallel=TRUE, is_negative_binomial=FALSE)
  stopCluster(cl)

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)

  # multiple chains in parallel using NB model
  cl <- makeCluster(2)
  registerDoParallel(cl)
  priors$overdispersion <- list(mean=5, sd=10)
  res <- mcmc(niterations=niter,
              snapshot_with_Rt_index_df,
              priors,
              serial_parameters,
              initial_cases_true,
              initial_reporting_parameters,
              initial_Rt,
              reporting_metropolis_parameters=list(mean_step=0.25, sd_step=0.1),
              serial_max=40, p_gamma_cutoff=0.99, maximise=FALSE,
              nchains=nchains, is_parallel=TRUE, is_negative_binomial=TRUE)
  stopCluster(cl)

  cases_df <- res$cases
  Rt_df <- res$Rt
  rep_df <- res$reporting
  overdispersion_df <- res$overdispersion

  expect_equal(max(cases_df$chain), nchains)
  expect_equal(max(Rt_df$chain), nchains)
  expect_equal(max(rep_df$chain), nchains)
  expect_equal(max(overdispersion_df$chain), nchains)
})

test_that("construct_w_matrix works ok", {
  w <- c(0.5, 0.25, 0.25)
  piece_width <- 5
  m_w <- construct_w_matrix(w, piece_width)

  wmax <- length(w)
  expect_equal(nrow(m_w), piece_width)
  expect_equal(ncol(m_w), wmax + piece_width - 1)

  row_sums <- rowSums(m_w)
  expect_true(sum(row_sums==1) == piece_width)
})

test_that("nb_log_likelihood_Rt_piece works as expected", {

  w <- c(0.5, 0.25, 0.25)
  Rt <- 2.2
  kappa <- 3.5
  onset_times <- c(4, 5)
  cases_df <- dplyr::tribble(
    ~time_onset, ~cases_true,
    1, 1,
    2, 3,
    3, 6,
    4, 5,
    5, 3
  )

  # check throws error if cases_df doesn't have onset times at bottom
  cases_df_wrong <- cases_df %>%
    dplyr::arrange(desc(time_onset))
  expect_error(nb_log_likelihood_Rt_piece(
    Rt, kappa, w, onset_times, cases_df_wrong)
  )

  # check works out logp pointwise
  logp <- nb_log_likelihood_Rt_piece(
    Rt, kappa, w, onset_times, cases_df
  )
  mu <- sum(Rt * w * rev(cases_df$cases_true[1:3]))
  logp_1 <- dnbinom(cases_df$cases_true[4], mu=mu, size=kappa,
                    log = TRUE)
  mu <- sum(Rt * w * rev(cases_df$cases_true[2:4]))
  logp_2 <- dnbinom(cases_df$cases_true[5], mu=mu, size=kappa,
                    log = TRUE)
  expect_equal(logp_1 + logp_2, logp)

  # check onset_time 1 contributes no logp
  onset_times <- c(1, 2)
  cases_df <- cases_df %>%
    dplyr::filter(time_onset %in% onset_times)
  logp <- nb_log_likelihood_Rt_piece(
    Rt, kappa, w, onset_times, cases_df
  )
  mu <- sum(Rt * w * c(cases_df$cases_true[1], 0, 0))
  logp_2 <- dnbinom(cases_df$cases_true[2], mu=mu, size=kappa,
                    log = TRUE)
  expect_equal(logp, logp_2)

})

test_that("sample_nb_Rt_piece works as expected", {

  prior_shape <- 1
  prior_rate <- 0.2

  w <- c(0.5, 0.25, 0.25)
  onset_times <- c(4, 5)
  cases_df <- dplyr::tribble(
    ~time_onset, ~cases_true,
    1, 10,
    2, 30,
    3, 60,
    4, 50,
    5, 30
  )

  posterior_shape <- prior_shape + sum(cases_df$cases_true[4:5])
  posterior_rate <- prior_rate + sum(w * rev(cases_df$cases_true[1:3])) +
    sum(w * rev(cases_df$cases_true[2:4]))

  ndraws <- 1000
  nresamples <- 10000

  # for large kappa the NB and Poisson distributions should be similar
  kappa <- 10000
  Rt_nb <- sample_nb_Rt_piece(
    prior_shape, prior_rate,
    posterior_shape, posterior_rate,
    kappa,
    w,
    onset_times,
    cases_df,
    ndraws,
    nresamples)
  expect_equal(length(Rt_nb), ndraws)
  Rt_poisson <- stats::rgamma(ndraws, posterior_shape, posterior_rate)
  expect_true(abs(mean(Rt_nb) - mean(Rt_poisson)) < 1)
  expect_true(abs(sd(Rt_nb) - sd(Rt_poisson)) < 0.1)

  # for small kappa the sd of NB should be greater
  kappa <- 1
  Rt_nb <- sample_nb_Rt_piece(
    prior_shape, prior_rate,
    posterior_shape, posterior_rate,
    kappa,
    w,
    onset_times,
    cases_df,
    ndraws,
    nresamples)
  expect_true(sd(Rt_nb) > sd(Rt_poisson))
})
