% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probabilities.R
\name{state_process_logp}
\alias{state_process_logp}
\title{Determines the probability of a given number of cases given a history of them}
\usage{
state_process_logp(cases_true, cases_history, Rt, serial_parameters)
}
\arguments{
\item{cases_true}{true case count originating on day_onset}

\item{cases_history}{a vector containing history of cases arranged from
recent to past}

\item{Rt}{effective reproduction number}

\item{serial_parameters}{named list of 'mean' and 'sd' of gamma distribution
characterising the serial interval distribution}
}
\value{
a log-probability or vector of log-probabilities
}
\description{
Calculates:
\deqn{Pois(cases_true_t|Rt * \sum_tau=1^t_max w_t cases_true_t-tau)}
}
