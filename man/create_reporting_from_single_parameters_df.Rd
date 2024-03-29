% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulated-data.R
\name{create_reporting_from_single_parameters_df}
\alias{create_reporting_from_single_parameters_df}
\title{Creates a reporting parameter tibble with a single set of reporting
parameters}
\usage{
create_reporting_from_single_parameters_df(time_onsets, reporting_parameters)
}
\arguments{
\item{time_onsets}{a vector of symptom onset times}

\item{reporting_parameters}{named list of 'mean' and 'sd' of gamma distribution
characterising the reporting delay distribution}
}
\value{
a tibble with three columns: 'time_onset', 'mean', 'sd'
}
\description{
Creates a reporting parameter tibble with a single set of reporting
parameters
}
