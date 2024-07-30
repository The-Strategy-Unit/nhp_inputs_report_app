








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

#' Create mixture distributions for each activity subset:
#'
#' Aggregates over the peers to create a mixture distribution for each activity
#' subset.
#'
#' @param data A dataframe of with the mu and sigma of the normal distribution
#' for each peer and activity subset.
#'
#' @return A list of mixture distributions for each activity subset.
get_mixture_distributions <- function(data, activity_subsets){

  mix_dists <- list()

  for (i in (activity_subsets)) {
    dist_list <- list()

    peers <- data |>
      dplyr::filter(activity_subset == i) |>
      dplyr::distinct(peer) |>
      dplyr::pull()

    for (j in (peers)) {
      norm_param <- data |>
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

  return(mix_dists)
}

#' Get percentiles from the mixture distributions.
#'
#' Gets the percentiles for the ECDF and PDF of each activity subset's mixture
#' distribution.
#'
#' @param data A list of mixture distributions for each activity subset.
#' @param activity_subsets A vector of the unique activity subsets.
#'
#' @return A long dataframe of the percentiles for the ECDF and PDF of each
#' activity subset's mixture distribution.
get_percentiles <- function(data, activity_subsets){

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
      ecdf_value = data[[i]]@p(q = seq(0, 100, 1))
    ) |>
      dplyr::mutate(pdf_value = ecdf_value - lag(ecdf_value, 1))

    peer_agg_ecdf_pdf <- peer_agg_ecdf_pdf |>
      dplyr::bind_rows(activity_subset_ecdf_pdf)

    rm(activity_subset_ecdf_pdf)

  }

  return(peer_agg_ecdf_pdf)

}

#' Get mu and sigma from aggregating over normal distributions.
#'
#' @param data  A dataframe with the mu and sigma of the normal distribution
#' for each peer and activity subset.
#'
#' @return A dataframe with mu, sigma and number of peers for each mixture
#' distribution.
#'
#' Note mean and sd of unweighted mixture distribution of Normal distributions
#' is:
#' mu_mix = sum(mu) / n
#' sd_mix = ( sum(mu^2 + sigma^2) / n ) ^ (1/2)
#' source : https://stats.stackexchange.com/questions/447626/mean-and-variance-of-a-mixture-distribution

get_mu_sigma <- function(data){

  summary <- data |>
    dplyr::summarise(
      mu = mean(mu),
      sd = (mean(mu ^ 2 + sigma ^ 2)) ^ (1 / 2),
      peers = dplyr::n(),
      .by = activity_subset
    )

  return(summary)
}

#' Get p10, p50 and p90 of mixture distributions.
#'
#' @param data A list of mixture distributions for each activity subset.
#' @param activity_subsets A vector of the unique activity subsets.
#'
#' @return A dataframe with p10, p50 and p90 of each mixture distribution.
get_p10_p50_p90 <- function(data, activity_subsets){

  peer_agg_p10_p50_p90 <- data.frame(
    activity_subset = character(),
    p10 = numeric(),
    p50 = numeric(),
    p90 = numeric()
  )

  for (i in (1:length(activity_subsets))) {
    activity_subset_p10_p50_p90 <- data.frame(
      activity_subset = activity_subsets[i],
      p10 = data[[i]]@q(p = 0.1),
      p50 = data[[i]]@q(p = 0.5),
      p90 = data[[i]]@q(p = 0.9)
    )

    peer_agg_p10_p50_p90 <- peer_agg_p10_p50_p90 |>
      dplyr::bind_rows(activity_subset_p10_p50_p90)

    rm(activity_subset_p10_p50_p90)

  }

  return(peer_agg_p10_p50_p90)
}

#' Summarise mixture distributions.
#'
#' Get mu, sigma, number of peers, p10, p50 and p90 for each mixture
#' distribution.
#'
#' @param normal_dists A dataframe with the mu and sigma of the normal
#' distribution for each peer and activity subset.
#' @param mix_dists A list of mixture distributions for each activity subset.
#' @param activity_subsets A vector of the unique activity subsets.
#'
#' @return A dataframe of distribution characteristics.
get_distribution_characteristics <- function(normal_dists,
                                             mix_dists,
                                             activity_subsets) {
  peer_agg_mu_sigma_n <- get_mu_sigma(normal_dists)

  peer_agg_p10_p50_p90 <- get_p10_p50_p90(mix_dists, activity_subsets)

  peer_agg_dist_summary <- peer_agg_mu_sigma_n |>
    dplyr::left_join(peer_agg_p10_p50_p90, dplyr::join_by(activity_subset))

  return(peer_agg_dist_summary)

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
activity_subsets <- normal_dists |>
  dplyr::distinct(activity_subset) |>
  dplyr::pull()

mix_dists <- get_mixture_distributions(normal_dists, activity_subsets)

# 4 Capture percentiles for ecdfs and pdfs ----
peer_agg_ecdf_pdf <- get_percentiles(mix_dists, activity_subsets)

# 5 Capture distribution characteristics ----
# mean, sd, n, p10, p50, p90
peer_agg_dist_summary <- get_distribution_characteristics(normal_dists,
                                                          mix_dists,
                                                          activity_subsets)

rm(activity_subsets)

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
