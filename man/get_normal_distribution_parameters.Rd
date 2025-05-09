% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_mixture_distributions.R
\name{get_normal_distribution_parameters}
\alias{get_normal_distribution_parameters}
\title{Get normal distribution parameters from lo and hi values.}
\usage{
get_normal_distribution_parameters(data)
}
\arguments{
\item{data}{A dataframe with lo and hi values where each row is a different
normal distribution.}
}
\value{
A dataframe with the mu and sigma of the normal distribution for each
 row.
}
\description{
Creates a dataframe of the normal distribution parameters for each pair of lo
(p10) and hi (p90) values given. Rows where point estimates (lo = hi) or
where default values (lo = 0 and hi = 1) are given are excluded.
}
