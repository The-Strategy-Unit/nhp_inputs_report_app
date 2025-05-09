# Functions to create and visualise mixture distributions for each strategy by
# aggregrating over all peers.

#' Summarise mixture distributions.
#'
#' Get mu, sigma, number of peers, p10, p50 and p90 for each mixture
#' distribution.
#'
#' @param normal_dists A dataframe with the mu and sigma of the normal
#' distribution for each peer and mitigator.
#' @param mix_dists A list of mixture distributions for each mitigator.
#' @param mitigators A vector of the unique mitigators.
#'
#' @return A dataframe of distribution characteristics.
get_distribution_characteristics <- function(normal_dists,
                                             mix_dists,
                                             mitigators) {
  peer_agg_mu_sigma_n <- get_mu_sigma(normal_dists)

  peer_agg_p10_p50_p90 <- get_p10_p50_p90(mix_dists, mitigators)

  peer_agg_dist_summary <- peer_agg_mu_sigma_n |>
    dplyr::left_join(peer_agg_p10_p50_p90, dplyr::join_by(mitigator_code))

  return(peer_agg_dist_summary)

}

#' Create mixture distributions for each mitigator:
#'
#' Aggregates over the peers to create a mixture distribution for each activity
#' subset.
#'
#' @param data A dataframe of with the mu and sigma of the normal distribution
#' for each peer and mitigator.
#'
#' @return A list of mixture distributions for each mitigator.
get_mixture_distributions <- function(data, mitigators){

  mix_dists <- list()

  for (i in (mitigators)) {
    dist_list <- list()

    peers <- data |>
      dplyr::filter(mitigator_code == i) |>
      dplyr::distinct(scheme_code) |>
      dplyr::pull()

    for (j in (peers)) {
      norm_param <- data |>
        dplyr::filter(mitigator_code == i, scheme_code == j) |>
        dplyr::select(mu, sigma)

      peer_dist <- distr::Norm(mean = norm_param$mu, sd = norm_param$sigma)

      dist_list <- append(dist_list, peer_dist)

      rm(peer_dist, norm_param)

    }

    mitigator_mix_dist <- distr::UnivarMixingDistribution(Dlist = dist_list)

    mix_dists <- append(mix_dists, mitigator_mix_dist)

    rm(mitigator_mix_dist, dist_list, peers)

  }

  return(mix_dists)
}


#' Get mu and sigma from aggregating over normal distributions.
#'
#' @param data  A dataframe with the mu and sigma of the normal distribution
#' for each peer and mitigator.
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
      .by = mitigator_code
    )

  return(summary)
}

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
    # ensure distinct combinations of mitigator and scheme
    dplyr::distinct(mitigator_code, scheme_code, .keep_all = TRUE) |>
    # convert hi and lo to numerator of percentage - e.g. 90 instead of 0.9
    dplyr::mutate(
      value_lo = value_lo * 100,
      value_hi = value_hi * 100
    ) |>
    dplyr::filter(
      !(value_lo == 0 & value_hi == 1 | # Exclude default values.
        value_lo == value_hi)# Exclude point estimates
    ) |>
    dplyr::mutate(
      mu = (value_lo + value_hi) / 2,
      sigma = (value_hi - mu) / stats::qnorm(p = 0.90, mean = 0, sd = 1)
    )

  return(normal_dists)

}

#' Get p10, p50 and p90 of mixture distributions.
#'
#' @param data A list of mixture distributions for each mitigator.
#' @param mitigators A vector of the unique mitigators.
#'
#' @return A dataframe with p10, p50 and p90 of each mixture distribution.
get_p10_p50_p90 <- function(data, mitigators){

  peer_agg_p10_p50_p90 <- data.frame(
    mitigator_code = character(),
    p10 = numeric(),
    p50 = numeric(),
    p90 = numeric()
  )

  for (i in (1:length(mitigators))) {
    mitigator_p10_p50_p90 <- data.frame(
      mitigator_code = mitigators[i],
      p10 = data[[i]]@q(p = 0.1),
      p50 = data[[i]]@q(p = 0.5),
      p90 = data[[i]]@q(p = 0.9)
    )

    peer_agg_p10_p50_p90 <- peer_agg_p10_p50_p90 |>
      dplyr::bind_rows(mitigator_p10_p50_p90)

    rm(mitigator_p10_p50_p90)

  }

  return(peer_agg_p10_p50_p90)
}

#' Get percentiles from the mixture distributions.
#'
#' Gets the percentiles for the ECDF and PDF of each mitigator's mixture
#' distribution.
#'
#' @param data A list of mixture distributions for each mitigator.
#' @param mitigators A vector of the unique mitigators.
#'
#' @return A long dataframe of the percentiles for the ECDF and PDF of each
#' mitigator's mixture distribution.
get_percentiles <- function(data, mitigators){

  peer_agg_ecdf_pdf <- data.frame(
    mitigator_code = character(),
    q = numeric(),
    ecdf_value = numeric(),
    pdf_value = numeric()
  )

  for (i in (1:length(mitigators))) {
    mitigator_ecdf_pdf <- data.frame(
      mitigator_code = mitigators[i],
      q = seq(0, 100, 1),
      ecdf_value = data[[i]]@p(q = seq(0, 100, 1))
    ) |>
      dplyr::mutate(pdf_value = ecdf_value - dplyr::lag(ecdf_value, 1))

    peer_agg_ecdf_pdf <- peer_agg_ecdf_pdf |>
      dplyr::bind_rows(mitigator_ecdf_pdf)

    rm(mitigator_ecdf_pdf)

  }

  return(peer_agg_ecdf_pdf)

}


#' Probability plot.
#'
#' @param data A dataframe where each row is a percentile of the ECDF and PDF of
#' a mixture distribution. Distribution characteristics and labels, groups are
#' also included.
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
    ggplot2::facet_wrap(ggplot2::vars(mitigator_label),
                        scale = 'free_y') +
    ggplot2::scale_y_continuous(name = y_axis_label,
                                limits = c(0, y_axis_max)) +
    ggplot2::scale_x_continuous(name = 'estimated percentage reduction') +
    ggplot2::labs(title = title,
                  caption = 'p10 and p90 (blue), mean (red)')

  return(plot)

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

#' Prepare data for plotting densities.
#'
#' Wrangles data to a dataframe for plotting ECDFs and PDFs
#'
#' @param peer_agg_ecdf_pdf A long dataframe of the percentiles for the ECDF and
#' PDF of each mitigator's mixture distribution.
#' @param peer_agg_dist_summary A dataframe of distribution characteristics.
#' @param strategy_lookup A dataframe of labels, types and groups for each
#' strategy.
#'
#' @return Tibble of data showing mixture distributions
wrangle_data_for_density_plots <- function(peer_agg_ecdf_pdf,
                                           peer_agg_dist_summary,
                                           strategy_lookup) {

  # make strategy lookup names tidy
  strategy_lookup <- strategy_lookup |> janitor::clean_names()

  data_wrangled <- peer_agg_ecdf_pdf |>
    dplyr::left_join(peer_agg_dist_summary, dplyr::join_by(mitigator_code)) |>
    dplyr::left_join(strategy_lookup, dplyr::join_by(mitigator_code)) |>
    dplyr::mutate(mitigator_label = paste(mitigator_name,
                                          ' (n = ', peers, ')',
                                          sep = ''))

  return(data_wrangled)

}


#' Get mixture distribution data
#'
#' Co-ordinates the production of mixture distribution data from the various
#' functions above.
#'
#' @param dat Tibble - of mitigator data as produced by `populate_table()` in `fct_tabulate.R`
#'
#' @return Tibble of mixture distribution data for each mitigator
#' @export
get_mixture_distributions_dat <- function(dat) {

  # prepare input data
  # ensure distinct combinations of mitigator and scheme to tackle an
  # issue with multiple values by the same scheme for the same mitigator
  dat <- dat |>
    dplyr::distinct(
      mitigator_code,
      scheme_code,
      .keep_all = TRUE
    )

  # get a vector of mitigator codes
  mitigators <- dat |>
    dplyr::distinct(mitigator_code) |>
    dplyr::pull()

  # estimate mean and sd for each mitigator
  normal_dists <- get_normal_distribution_parameters(data = dat)

  # extract mitigator details from dat
  strategy_lookup <- dat |>
    dplyr::select(dplyr::starts_with('mitigator')) |>
    dplyr::distinct()

  # get mixture distributions for each mitigator
  mix_dists <- get_mixture_distributions(
    data = normal_dists,
    mitigators = mitigators
  )

  # capture percentiles for ecdfs and pdfs
  peer_agg_ecdf_pdf <- get_percentiles(
    data = mix_dists,
    mitigators = mitigators
  )

  # capture distribution characteristics
  peer_agg_dist_summary <- get_distribution_characteristics(
    normal_dists = normal_dists,
    mix_dists = mix_dists,
    mitigators = mitigators
  )

  # get data for plotting
  data_for_plotting <- wrangle_data_for_density_plots(
    peer_agg_ecdf_pdf = peer_agg_ecdf_pdf,
    peer_agg_dist_summary = peer_agg_dist_summary,
    strategy_lookup = strategy_lookup
  )
}
