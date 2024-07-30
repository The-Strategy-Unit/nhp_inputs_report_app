








# Functions ----

#' Get normal distribution parameters from lo and hi values.
#'
#' Creates a dataframe of the normal distribution parameters for each pair of lo
#' (p10) and hi (p90) values given. Rows where point estimates (lo = hi) or
#' where default values (lo = 0 and hi = 1) are given are excluded.
#'
#' @param data A dataframe with lo and hi values where each row is a different
#' normal distribution.
#'
#' @return A dataframe with the mu and sigma of the normal distribution for each
#'  row.
get_normal_distribution_parameters <- function(data) {
  normal_dists <- data |>

    dplyr::filter(!(lo == 0 & hi == 1 | # Exclude default values.
                      lo == hi)# Exclude point estimates
                  ) |>
    dplyr::mutate(mu = (lo + hi) / 2,
                  sigma = (hi - mu) / qnorm(p = 0.90, mean = 0, sd = 1))

  return(normal_dists)

}



#' Modifies theme of ECDF and PDF plot.
#'
#' @param plot A plot of an ECDF or PDF.
#' @param type Either `"ecdf"` or `"pdf"` to get the modifications for the
#' empirical cumulative distribution functions or probability density functions,
#' respectively.
#'
#' @return A plot with theme modifiers based on whether the plot was is an ECDF
#' or PDF.
modify_theme <- function(plot, type) {
  if (type == "ecdf") {
    plot <- plot
  } else {
    plot <- plot +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }

  return(plot)
}

#' Probability plot.
#'
#' @param data A dataframe where each row is ???????????????????????????????????????????????
#' @param type Either `"ecdf"` or `"pdf"` to get the empirical cumulative
#' distribution functions or probability density functions, respectively.
#'
#' @return A plot of the ECDF or PDF.
get_probability_plot <- function(data, type) {
  if (type == "ecdf") {
    title <- "Empirical cumulative distribution functions | aggregated peer opinions"
    y_axis_label <- "probability"
    y_axis_max <- 1
  } else {
    title <- "Probability density functions | aggregated peer opinions"
    y_axis_label <- "probability density"
    y_axis_max <- NA_real_
  }

  y <- paste0(type, "_value")

  plot <- data |>
    ggplot2::ggplot() |>
    modify_theme(type) +
    ggplot2::geom_line(ggplot2::aes(x = q, y = !!rlang::sym(y))) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = p10), colour = 'blue') +
    ggplot2::geom_vline(ggplot2::aes(xintercept = p90), colour = 'blue') +
    ggplot2::geom_vline(ggplot2::aes(xintercept = mu), colour = 'red') +
    ggplot2::facet_wrap(ggplot2::vars(activity_subset_label),
                        scale = 'free_y') +
    ggplot2::scale_y_continuous(name = y_axis_label,
                                limits = c(0, y_axis_max)) +
    ggplot2::scale_x_continuous(name = 'estimated percentage reduction') +
    ggplot2::labs(title = title,
                  subtitle = 'inpatient admissions avoidance',
                  caption = 'p10 and p90 (blue), mean (red)')

  return(plot)

}



# 0 Setup ----
library(distr)
library(tidyverse)

# 1 Read data ----
# Start with a dataframe that has at least these columns:
# peer - fct - usually organisation code
# strategy - chr - the coder-friendly names of the mitigators
# lo - dbl - p10 of the distribution
# hi - dbl - p90 of the distribution
# note lo and hi should be the numerator of the percentage - for example 90
# instead of 0.9

data <- readRDS(file = "input_data.rds") |>
  rename(activity_subset = strategy)

strategy_lookup <- read.csv('strategy_lookup.csv', header = TRUE)

# 2 Wrangle data ----
normal_dists <- get_normal_distribution_parameters(data)

# 3 Aggregate estimates ----
# Create mixture distributions for each activity subset:
activity_subsets <- normal_dists |>
  dplyr::distinct(activity_subset) |>
  dplyr::pull()

mix_dists <- list()

for (i in (activity_subsets)) {
  dist_list <- list()

  peers <- normal_dists |>
    dplyr::filter(activity_subset == i) |>
    dplyr::distinct(peer) |>
    dplyr::pull()

  for (j in (peers)) {
    norm_param <- normal_dists |>
      dplyr::filter(activity_subset == i, peer == j) |>
      dplyr::select(mu, sigma)


    peer_dist <- distr::Norm(mean = norm_param$mu, sd = norm_param$sigma)

    dist_list <- append(dist_list, peer_dist)

    rm(peer_dist, norm_param)

  }

  activity_subset_mix_dist <- distr::UnivarMixingDistribution(Dlist = dist_list)

  mix_dists <- append(mix_dists, activity_subset_mix_dist)

  rm(activity_subset_mix_dist, dist_list, peers)

}

# 4 Capture percentiles for ecdfs and pdfs ----
peer_agg_ecdf_pdf <- data.frame(
  activity_subset = character(),
  q = numeric(),
  ecdf_value = numeric(),
  pdf_value = numeric()
)


for (i in (1:length(activity_subsets))) {
  activity_subset_ecdf_pdf <- data.frame(
    activity_subset = activity_subsets[i],
    q = seq(0, 100, 1),
    ecdf_value = mix_dists[[i]]@p(q = seq(0, 100, 1))
  ) |>
    dplyr::mutate(pdf_value = ecdf_value - lag(ecdf_value, 1))

  peer_agg_ecdf_pdf <- peer_agg_ecdf_pdf |>
    dplyr::bind_rows(activity_subset_ecdf_pdf)

  rm(activity_subset_ecdf_pdf)

}

# 5 Capture distribution characteristics ----
# mean, sd, p10, p50, p90

# note mean and sd of unweighted mixture distribution of Normal distributions is
# mu_mix = sum(mu) / n
# sd_mix = ( sum(mu^2 + sigma^2) / n ) ^ (1/2)
# source : https://stats.stackexchange.com/questions/447626/mean-and-variance-of-a-mixture-distribution

peer_agg_mu_sigma_n <- normal_dists |>
  dplyr::summarise(
    mu = mean(mu),
    sd = (mean(mu ^ 2 + sigma ^ 2)) ^ (1 / 2),
    peers = dplyr::n(),
    .by = activity_subset
  )

peer_agg_p10_p50_p90 <- data.frame(
  activity_subset = character(),
  p10 = numeric(),
  p50 = numeric(),
  p90 = numeric()
)

for (i in (1:length(activity_subsets))) {
  activity_subset_p10_p50_p90 <- data.frame(
    activity_subset = activity_subsets[i],
    p10 = mix_dists[[i]]@q(p = 0.1),
    p50 = mix_dists[[i]]@q(p = 0.5),
    p90 = mix_dists[[i]]@q(p = 0.9)
  )

  peer_agg_p10_p50_p90 <- peer_agg_p10_p50_p90 |>
    dplyr::bind_rows(activity_subset_p10_p50_p90)

  rm(activity_subset_p10_p50_p90)

}

peer_agg_dist_summary <- peer_agg_mu_sigma_n |>
  dplyr::left_join(peer_agg_p10_p50_p90, dplyr::join_by(activity_subset))

rm(peer_agg_mu_sigma_n, peer_agg_p10_p50_p90, activity_subsets)

# Save data:
#saveRDS(peer_agg_dist_summary, file = "mixture_distributions_output.rds")

# 6 Visualise results ----
# ecds and pdfs for admission avoidance

data_admission_avoidance <- peer_agg_ecdf_pdf |>
  dplyr::left_join(peer_agg_dist_summary, dplyr::join_by(activity_subset)) |>
  dplyr::left_join(strategy_lookup, dplyr::join_by(activity_subset == strategy)) |>
  dplyr::mutate(activity_subset_label = paste(strategyLabel, ' (n = ', peers, ')', sep = '')) |>
  dplyr::filter(strategyType == 'inpatient admission avoidance')



# Note I've put a filter on strategyGroup so the plots can be seen more clearly
# in this toy example, but this can be removed / amended as needed:
data_admission_avoidance |>
  dplyr::filter(strategyGroup ==
                  "Hospital activity amenable to public health interventions") |>
  get_probability_plot(type = "ecdf")

data_admission_avoidance |>
  dplyr::filter(strategyGroup ==
                  "Hospital activity amenable to public health interventions") |>
  get_probability_plot(type = "pdf")
