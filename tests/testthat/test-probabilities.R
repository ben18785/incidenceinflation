test_that("observation_process_single_logp returns reasonable log-probs", {
  reporting_parameters <- list(mean=5, sd=4)
  day_2 <- 4
  day_1 <- 2
  day_onset <- 0
  I_true <- 10
  I_day_1 <- 3

  # I_day_2 > I_true should be impossible
  I_day_2 <- 11
  logp <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  expect_equal(-Inf, logp)

  # I_day_2 < I_day_1 should be impossible
  I_day_2 <- 2
  logp <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  expect_equal(-Inf, logp)

  # I_day_2 near I_true should be more likely with shorter delays
  I_day_2 <- 10
  logp <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  reporting_parameters1 <- list(mean=1, sd=1)
  logp1 <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                    day_2, day_1, day_onset, reporting_parameters1)
  expect_true(logp1 > logp)

  # prob 1 if I_day_1=I_true
  I_day_2 <- 10
  I_day_1 <- 10
  logp <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                   day_2, day_1, day_onset, reporting_parameters)
  expect_equal(0, logp)
})

test_that("observation_process_single_logp can return vectors", {
  reporting_parameters <- list(mean=5, sd=4)
  day_2 <- 4
  day_1 <- 2
  day_onset <- 0
  I_true <- 10:20
  I_day_1 <- 3
  I_day_2 <- 5
  logps <- observation_process_single_logp(I_true, I_day_2, I_day_1,
                                          day_2, day_1, day_onset, reporting_parameters)
  n <- length(I_true)
  expect_equal(length(logps), n)
  log_p_long_way <- vector(length=n)
  for(i in 1:n)
    log_p_long_way[i] <- observation_process_single_logp(
      I_true[i], I_day_2, I_day_1, day_2, day_1, day_onset, reporting_parameters)
  expect_true(all.equal(logps, log_p_long_way))
})

test_that("observation_process_logp returns reasonable log prob values", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 20, 21),
                                      cases_reported=c(3, 5, 7))
  reporting_parameters <- list(mean=1, sd=1)
  logp <- observation_process_logp(observation_df=observation_matrix,
                                   cases_true = 10,
                                   day_onset = 1,
                                   reporting_parameters=reporting_parameters)
  expect_true(logp < 0)

  # test that prob declines as delay increases
  observation_matrix <- dplyr::tibble(time_reported=c(1, 20, 21),
                                      cases_reported=c(1, 1, 1))
  logp1 <- observation_process_logp(observation_df=observation_matrix,
                                   cases_true = 10,
                                   day_onset = 0,
                                   reporting_parameters=reporting_parameters)
  expect_true(logp1 < logp)

  observation_matrix <- dplyr::tibble(time_reported=c(1, 1, 5),
                                      cases_reported=c(3, 5, 7))
  expect_error(observation_process_logp(observation_df=observation_matrix,
                                   cases_true = 10,
                                   day_onset = 1,
                                   reporting_parameters=reporting_parameters))
  observation_matrix <- dplyr::tibble(time_reported=c(1, 2, 5),
                                      cases_reported=c(3, 2, 7))
  expect_error(observation_process_logp(observation_df=observation_matrix,
                                        cases_true = 10,
                                        day_onset = 1,
                                        reporting_parameters=reporting_parameters))
})

test_that("state_process_logp produces reasonable log prob values", {
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  w <- weights_series(t_max, s_params)
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 1
  expect_equal(state_process_logp(cases_true, cases_history, Rt, s_params), -1)

  Rt <- 2
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 2
  expect_equal(state_process_logp(cases_true, cases_history, Rt, s_params), log(2 * exp(-2)))
})

test_that("state_process_logp produces works fine with vectors of true cases", {
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  w <- weights_series(t_max, s_params)
  a_sum <- Rt * sum(w * cases_history)
  cases_true <- 1:5
  logps <- state_process_logp(cases_true, cases_history, Rt, s_params)
  n <- length(cases_true)
  log_p_long_way <- vector(length=n)
  for(i in 1:n)
    log_p_long_way[i] <- state_process_logp(cases_true[i], cases_history, Rt, s_params)
  expect_true(all.equal(logps, log_p_long_way))
})

test_that("conditional_cases_logp returns reasonable values", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=1, sd=1)
  cases_true <- 80
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  day_onset <- 0
  logp <- conditional_cases_logp(observation_df=observation_matrix,
                         cases_true=cases_true,
                         cases_history=cases_history,
                         Rt=Rt,
                         day_onset=day_onset,
                         serial_parameters=s_params,
                         reporting_parameters=reporting_parameters)
  expect_true(logp < 0)
  reporting_parameters <- list(mean=20, sd=10)
  logp1 <- conditional_cases_logp(observation_df=observation_matrix,
                                 cases_true=cases_true,
                                 cases_history=cases_history,
                                 Rt=Rt,
                                 day_onset=day_onset,
                                 serial_parameters=s_params,
                                 reporting_parameters=reporting_parameters)
  expect_true(logp1 > logp)
})

test_that("conditional_cases_logp works fine with vectorised true cases", {
  observation_matrix <- dplyr::tibble(time_reported=c(1, 3, 5),
                                      cases_reported=c(1, 1, 1))
  reporting_parameters <- list(mean=1, sd=1)
  cases_true <- 10:20
  s_params <- list(mean=5, sd=2)
  Rt <- 1
  t_max <- 20
  cases_history <- rep(1, t_max)
  day_onset <- 0
  logps <- conditional_cases_logp(observation_df=observation_matrix,
                                 cases_true=cases_true,
                                 cases_history=cases_history,
                                 Rt=Rt,
                                 day_onset=day_onset,
                                 serial_parameters=s_params,
                                 reporting_parameters=reporting_parameters)
  n <- length(cases_true)
  log_p_long_way <- vector(length=n)
  for(i in 1:n)
    log_p_long_way[i] <- conditional_cases_logp(observation_df=observation_matrix,
                                                cases_true=cases_true[i],
                                                cases_history=cases_history,
                                                Rt=Rt,
                                                day_onset=day_onset,
                                                serial_parameters=s_params,
                                                reporting_parameters=reporting_parameters)
  expect_true(all.equal(logps, log_p_long_way))
})

days_total <- 100
Rt_1 <- 1.5
Rt_2 <- 1.0
Rt_3 <- 1.3
v_Rt <- c(rep(Rt_1, 40), rep(Rt_2, 20), rep(Rt_3, 40))
Rt_function <- stats::approxfun(1:days_total, v_Rt)
s_params <- list(mean=5, sd=3)
r_params <- list(mean=10, sd=3)
kappa <- 1000
days_total <- 3
df <- generate_snapshots(days_total, Rt_function, s_params, r_params,
                         kappa=kappa)

test_that("observation_process_all_times_logp works ok", {
  logp <- observation_process_all_times_logp(df, r_params)
  df1 <- df %>% dplyr::filter(time_onset==1)
  logp1 <- observation_process_logp(df1, df1$cases_true[1],
                                    1, r_params)
  df1 <- df %>% dplyr::filter(time_onset==2)
  logp2 <- observation_process_logp(df1, df1$cases_true[1],
                                    2,r_params)
  expect_equal(logp, logp1 + logp2)

  # check -inf if negative mu or sd
  r_params <- list(mean = -1 , sd = 2)
  logp <- observation_process_all_times_logp(df, r_params)
  expect_equal(logp, -Inf)
  r_params <- list(mean = 1 , sd = -2)
  logp <- observation_process_all_times_logp(df, r_params)
  expect_equal(logp, -Inf)
  df <- dplyr::tribble(
    ~time_onset, ~time_reported, ~cases_reported, ~cases_true,
    1, 1, 0, 49,
    1, 2, 0, 49,
    1, 3, 0, 50,
    2, 2, 0, 55,
    2, 3, 0, 55,
    3, 3, 0, 54
  )
  r_params <- list(mean = 1 , sd = 2)
  expect_error(
    observation_process_all_times_logp(df, r_params)
  )
})
