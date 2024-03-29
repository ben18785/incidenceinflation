% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{sample_nb_Rt_piece}
\alias{sample_nb_Rt_piece}
\title{Uses importance resampling to infer a posterior over Rt under a negative
binomial renewal model}
\usage{
sample_nb_Rt_piece(
  prior_shape,
  prior_rate,
  posterior_shape,
  posterior_rate,
  kappa,
  w,
  onset_times,
  cases_df,
  ndraws,
  nresamples
)
}
\arguments{
\item{prior_shape}{Rt prior shape parameter}

\item{prior_rate}{Rt prior rate parameter}

\item{posterior_shape}{Rt posterior shape parameter}

\item{posterior_rate}{Rt posterior rate parameter}

\item{kappa}{overdispersion parameter}

\item{w}{weights corresponding to the generation times}

\item{onset_times}{the onset times corresponding to the piece}

\item{cases_df}{a tibble with 'cases_true' as a column which has been ordered
so that the latest onset times are at the bottom}

\item{ndraws}{number of draws of Rt to return}

\item{nresamples}{number of resamples used to calculate weights for}
}
\value{
a vector of Rt draws
}
\description{
Uses importance resampling to infer a posterior over Rt under a negative
binomial renewal model
}
