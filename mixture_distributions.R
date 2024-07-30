# Script to create and visualise mixture distributions for each mitigator using
# example data.

# 0 Setup ----
library(distr)
library(tidyverse)

source("mixture_distributions_functions.R")

# 1 Read data ----
# Start with a dataframe that has at least these columns:
# peer - fct - usually organisation code
# strategy - chr - the coder-friendly names of the mitigators
# lo - dbl - p10 of the distribution
# hi - dbl - p90 of the distribution
# note lo and hi should be the numerator of the percentage - for example 90
# instead of 0.9

data <- readRDS(file = "input_data.rds") |>
  rename(mitigator = strategy)

strategy_lookup <- read.csv('strategy_lookup.csv', header = TRUE)

# 2 Wrangle data ----
normal_dists <- get_normal_distribution_parameters(data)

# 3 Aggregate estimates ----
mitigators <- normal_dists |>
  dplyr::distinct(mitigator) |>
  dplyr::pull()

mix_dists <- get_mixture_distributions(normal_dists, mitigators)

# 4 Capture percentiles for ecdfs and pdfs ----
peer_agg_ecdf_pdf <- get_percentiles(mix_dists, mitigators)

# 5 Capture distribution characteristics ----
peer_agg_dist_summary <- get_distribution_characteristics(normal_dists,
                                                          mix_dists,
                                                          mitigators)

# Save data:
#saveRDS(peer_agg_dist_summary, file = "mixture_distributions_output.rds")

# 6 Visualise results ----
data_for_plotting <- wrangle_data_for_density_plots(peer_agg_ecdf_pdf,
                                                    peer_agg_dist_summary,
                                                    strategy_lookup)

# Note there's a filter on strategyGroup so the plots can be seen more clearly
# in this toy example, but this can be removed / amended as needed:
data_for_plotting |>
  dplyr::filter(strategyGroup ==
                "Hospital activity amenable to public health interventions") |>
  get_probability_plot(type = "ecdf")

data_for_plotting |>
  dplyr::filter(strategyGroup ==
                "Hospital activity amenable to public health interventions") |>
  get_probability_plot(type = "pdf")
