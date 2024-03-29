% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probabilities.R
\name{observation_process_single_logp}
\alias{observation_process_single_logp}
\title{Calculates observation process log probability density for a pair of days}
\usage{
observation_process_single_logp(
  cases_true,
  cases_day_2,
  cases_day_1,
  day_2,
  day_1,
  day_onset,
  reporting_parameters
)
}
\arguments{
\item{cases_true}{true case count(s) originating on day_onset: can be a single value or vector}

\item{cases_day_2}{number of cases arising on day_onset observed by day_2}

\item{cases_day_1}{number of cases arising on day_onset observed by day_1 < day_2}

\item{day_2}{day of last observation: day_2 > day_1 >= day_onset}

\item{day_1}{day of first observation: day_1 >= day_onset}

\item{day_onset}{day when case originates}

\item{reporting_parameters}{named list of 'mean' and 'sd' of gamma distribution
characterising the reporting delay distribution}
}
\value{
a log-probability or vector of log-probabilities
}
\description{
The probability of observations is a function of the number of cases remaining
to be reported:
\deqn{I_remaining = I_true - I_day_1}
and the number of cases observed between day_1 and day_2:
\deqn{I_obs = I_day_2 - I_day_1}
It is given by:
\deqn{binomial(I_obs|I_remaining, p_detect)}
where p_detect is the probability a thus undetected case is detected between
day_1 and day_2.
}
