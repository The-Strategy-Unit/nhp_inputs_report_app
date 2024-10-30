#' Generate a Simple Demo Dataset
#' @param n_schemes Integer.
#' @param n_mitigators Integer.
#' @param n_mitigator_subsets Integer.
#' @param prop_selected Numeric.
#' @param seed Integer.
#' @return A data.frame.
#' @noRd
generate_test_dataset <- function(
    n_schemes = 20,
    n_mitigators = 50,
    n_mitigator_subsets = 8,
    prop_selected = 0.6,
    seed = 123
) {

  scheme <- paste(
    "Scheme",
    stringr::str_pad(1:n_schemes, width = 2, side = "left", pad = "0")
  )

  mitigator <- paste(
    "Mitigator",
    stringr::str_pad(1:n_mitigators, width = 2, side = "left", pad = "0")
  )

  mitigator_lookup <- mitigator |>
    split(
      cut(
        seq_along(mitigator),
        n_mitigator_subsets,
        labels = paste("Group", 1:n_mitigator_subsets)
      )
    ) |>
    tibble::enframe(name = "mitigator_subset", value = "mitigator") |>
    tidyr::unnest(mitigator)

  combos <- tidyr::crossing(scheme, mitigator) |>
    dplyr::left_join(mitigator_lookup, by = "mitigator")

  nrows <- nrow(combos)

  withr::with_seed(
    seed = seed,
    code = {
      combos |>
        dplyr::mutate(
          lo = runif(nrows),
          hi = pmin(lo + runif(nrows, 0.25, 0.5), 1),
          mid = (lo + hi) / 2,
          years = sample(20:30, nrows, replace = TRUE)
        ) |>
        dplyr::slice_sample(prop = prop_selected)
    }
  )

}

#' Generate a Lookup of Schemes to Peer set
#' @param test_dataset A data.frame
#' @param n_peers Integer.
#' @return A data.frame.
#' @noRd
generate_peer_set <- function(test_dataset, n_peers = 5) {

  schemes <- test_dataset |>
    dplyr::distinct(scheme) |>
    dplyr::arrange(scheme) |>
    dplyr::pull(scheme)

  peer_set <- vector("list", length = length(schemes)) |>
    setNames(schemes)

  for (i in seq_along(schemes)) {
    schemes_to_sample <- schemes[schemes != schemes[i]]
    peer_set[[i]] <- sample(schemes_to_sample, n_peers, replace = FALSE)
  }

  peer_set |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(
      tidyselect::everything(),
      names_to = "scheme",
      values_to = "peer"
    ) |>
    dplyr::arrange(scheme, peer)

}
