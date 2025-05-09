% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_plots.R
\name{plot_trendline_comparison}
\alias{plot_trendline_comparison}
\title{Plot an individual trendline comparison plot}
\usage{
plot_trendline_comparison(
  plot_data,
  dat_lu,
  mitigator_codes,
  focal_scheme_code,
  show_other_schemes = TRUE,
  show_horizon_timeline = TRUE,
  show_horizon_overlay = TRUE,
  show_prebaseline_average = TRUE,
  facet_columns = 1,
  facet_height_px = 400,
  facet_count = 1,
  x_axis_min = 2010,
  return_data = FALSE
)
}
\arguments{
\item{plot_data}{Tibble - historical rate of activity for the mitigator, as produced within `plot_facetted_trendlines()` in `fct_plots.R`}

\item{dat_lu}{Tibble - lookup information for schemes and mitigators, as produced within `plot_facetted_trendlines()` in `fct_plots.R`}

\item{mitigator_codes}{Character vector - the `mitigator_code` to produce the plot for}

\item{focal_scheme_code}{Character vector - the focal `scheme_code` to produce the plot for}

\item{show_other_schemes}{Boolean (default = TRUE) plot time series for other schemes}

\item{show_horizon_timeline}{Boolean (default = TRUE) plot the predicted activity on the timeline at the horizon year}

\item{show_horizon_overlay}{Boolean (default = TRUE) plot the predicted activity as an overlay over the historical time series plots}

\item{show_prebaseline_average}{Boolean (default = TRUE) show the pre-baseline average (mean) with a range of two standard deviations above and below}

\item{facet_columns}{Integer (default = 1) the number of columns to facet the baseline plot}

\item{facet_height_px}{Integer (default = 400) the height of the plot}

\item{facet_count}{Integer (default = 1) the number of individual plots to be combined}

\item{x_axis_min}{Integer (default = 2010) the minimum year to show on the x-axis to help coordinate the x-axes for each plot}

\item{return_data}{Boolean (default = FALSE) TRUE = return a list object of tibbles used in the production of the plot - FOR TROUBLESHOOTING PURPOSES ONLY}
}
\value{
{plotly} plot (or list of tibbles if return_data == TRUE)
}
\description{
Produces a single {plotly} object showing the trendline for the specified
scheme and mitigator.
}
\details{
Note - this function is designed for use with `plot_facetted_trendlines()`
which coordinates the production of plots for multiple mitigators.
The use of the `subplot()` features of {plotly} mean that the overall plot
height is specified in each subplot's definition.
}
