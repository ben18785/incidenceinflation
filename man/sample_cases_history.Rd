% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{sample_cases_history}
\alias{sample_cases_history}
\title{Draws a possible history (or histories) of cases}
\usage{
sample_cases_history(
  observation_onset_df,
  max_cases,
  Rt_function,
  serial_parameters,
  reporting_parameters,
  p_gamma_cutoff = 0.99,
  maximise = FALSE
)
}
\arguments{
\item{observation_onset_df}{a tibble with three columns: time_onset, time_reported, cases_reported}

\item{max_cases}{maximum possible cases thought to arise on a given day}

\item{Rt_function}{takes day as an input and outputs an Rt value}

\item{serial_parameters}{named list of 'mean' and 'sd' of gamma distribution
characterising the serial interval distribution}

\item{reporting_parameters}{named list of 'mean' and 'sd' of gamma distribution
characterising the reporting delay distribution}

\item{p_gamma_cutoff}{a p value (0 <= p <= 1) indicating the threshold above which
we deem certainty}

\item{maximise}{rather than sample a case count give the case count with the
maximum probability (by default is FALSE)}
}
\value{
a tibble with an extra cases_estimated column
}
\description{
The distribution being drawn from at each time t is given by:
\deqn{p(cases_true_t|data, Rt, reporting_params, serial_params) \propto p(data|cases_true, reporting_params)
 p(cases_true_t|cases_true_t_1, cases_true_t_2, ..., Rt, serial_params)}
}
