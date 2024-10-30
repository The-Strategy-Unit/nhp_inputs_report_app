extract_params <- function(params, runs_meta) {

  possibly_report_params_table <- purrr::possibly(report_params_table)

  activity_avoidance <- params |>
    purrr::map(possibly_report_params_table, "activity_avoidance") |>
    purrr::list_rbind()

  efficiencies <- params |>
    purrr::map(possibly_report_params_table, "efficiencies") |>
    purrr::list_rbind()

  runs_meta <- runs_meta |> dplyr::select(dataset, scenario, run_stage)

  activity_avoidance |>
    dplyr::bind_rows(efficiencies) |>
    dplyr::mutate(
      peer_year = paste0(
        peer,
        "_", stringr::str_sub(baseline_year, 3, 4),
        "_", stringr::str_sub(horizon_year, 3, 4)
      )
    ) |>
    dplyr::left_join(runs_meta, by = dplyr::join_by("peer" == "dataset")) |>
    correct_day_procedures()

}

correct_day_procedures <- function(x) {

  # Identify pairs of bads/day_procedures mitigators with flag
  flagged <- x |>
    dplyr::mutate(
      mitigator_code_flag = dplyr::case_when(
        stringr::str_detect(
          strategy,
          "^bads_daycase$|^day_procedures_usually_dc$"  # old name/new name
        ) ~ "IP-EF-005",  # might as well flag with the mitigator code
        stringr::str_detect(
          strategy,
          "^bads_daycase_occasional$|^day_procedures_occasionally_dc$"
        ) ~ "IP-EF-006",
        stringr::str_detect(
          strategy,
          "^bads_outpatients$|^day_procedures_usually_op$"
        ) ~ "IP-EF-007",
        stringr::str_detect(
          strategy,
          "^bads_outpatients_or_daycase$|^day_procedures_occasionally_op$"
        ) ~ "IP-EF-008",
        .default = NA_character_
      )
    )

  # Identify where a peer has more than one instance of the code, i.e. the
  # mitigator is represented by both a bads and a day_procedures version. We'll
  # use this info to filter out the bads version.
  dupes <- flagged |>
    dplyr::count(peer, mitigator_code_flag) |>
    tidyr::drop_na(mitigator_code_flag) |>
    dplyr::filter(n > 1)

  # Remove bads mitigators if there's a day_procedures replacement for it
  for (i in seq(nrow(dupes))) {
    flagged <- flagged |>
      dplyr::filter(
        !(peer == dupes[[i, "peer"]] &
            mitigator_code_flag == dupes[[i, "mitigator_code_flag"]] &
            stringr::str_detect(strategy, "^bads_"))
      )
  }

  # Remaining bads mitigators clearly don't have a replacement day_procedures
  # version so we can just rename these ones.
  flagged |>
    dplyr::mutate(
      strategy = dplyr::case_match(
        strategy,
        "bads_daycase" ~ "day_procedures_usually_dc",
        "bads_daycase_occasional" ~ "day_procedures_occasionally_dc",
        "bads_outpatients" ~ "day_procedures_usually_op",
        "bads_outpatients_or_daycase" ~ "day_procedures_occasionally_op",
        .default = strategy
      )
    ) |>
    dplyr::select(-mitigator_code_flag)  # remove helper column

}

report_params_table <- function(
    p,  # a single scheme's params
    parameter = c("activity_avoidance", "efficiencies")
) {

  parameter_data <- p[[parameter]]

  time_profiles <- p[["time_profile_mappings"]][[parameter]] |>
    purrr::map(unlist) |>
    purrr::map(tibble::enframe, "strategy", "time_profile") |>
    data.table::rbindlist(idcol = "activity_type") |>
    dplyr::tibble()

  parameter_data |>
    purrr::map_depth(2, "interval") |>
    purrr::map(tibble::enframe, "strategy") |>
    dplyr::bind_rows(.id = "activity_type") |>
    tidyr::unnest_wider("value", names_sep = "_") |>
    dplyr::left_join(
      time_profiles,
      by = dplyr::join_by("activity_type", "strategy")
    ) |>
    dplyr::arrange("activity_type_name", "mitigator_name") |>
    dplyr::mutate(
      parameter = parameter,
      peer = p[["dataset"]],
      baseline_year = p[["start_year"]],
      horizon_year = p[["end_year"]]
    )

}

prepare_skeleton_table <- function(extracted_params) {

  strategies <- extracted_params |>
    dplyr::distinct(activity_type, strategy, parameter)

  tidyr::expand_grid(
    "strategy" = unique(extracted_params[["strategy"]]),
    "peer_year" = unique(extracted_params[["peer_year"]])
  ) |>
    dplyr::left_join(strategies, by = dplyr::join_by(strategy)) |>
    dplyr::select(peer_year, activity_type, parameter, strategy) |>
    dplyr::arrange(peer_year, activity_type, parameter, strategy)

}

#' Populate dat table
#'
#' Creates a single fact table from `skeleton_table` by left-joining with
#' `extracted_params`, `trust_code_lookup`, `mitigator_lookup` and
#' `nee_results`.
#'
#' @param skeleton_table Tibble - output of `prepare_skeleton_table()`
#' @param extracted_params Tibble - output of `extract_params()`
#' @param trust_code_lookup Tibble - from Azure file `nhp-scheme-lookup.csv`
#' @param mitigator_lookup Tibble - from Azure file `mitigator-lookup.csv`
#' @param nee_results Tibble - from fct_read.R of `nee_table.rds`
#'
#' @return Tibble of data
#' @export
populate_table <- function(
    skeleton_table,
    extracted_params,
    trust_code_lookup,
    mitigator_lookup,
    nee_results
) {

  data_joined <- skeleton_table |>
    dplyr::left_join(
      extracted_params,
      by = dplyr::join_by(peer_year, activity_type, parameter, strategy)
    ) |>
    dplyr::left_join(
      trust_code_lookup |>
        dplyr::mutate(
          scheme_code = `Trust ODS Code`,
          scheme_name = `Name of Hospital site`,
          .keep = "none"
        ),
      by = dplyr::join_by(peer == scheme_code)
    ) |>
    dplyr::left_join(
      mitigator_lookup,
      by = dplyr::join_by(strategy == "Strategy variable")
    ) |>
    dplyr::left_join(
      nee_results,
      by = dplyr::join_by(strategy == param_name)
    )

  data_prepared <- data_joined |>
    dplyr::filter(
      (value_1 <= 1 & value_2 <= 1) |
        is.na(value_1) & is.na(value_2)
    ) |>
    dplyr::mutate(
      .keep = "none",
      # schemes
      scheme_name = dplyr::case_when(
        # identify schemes whose mitigators are not yet finalised
        !is.na(scheme_name) & stringr::str_starts(
          string = run_stage |> stringr::str_to_lower(),
          pattern = 'final',
          negate = TRUE
        ) ~ glue::glue('{scheme_name} ✏️'),
        .default = scheme_name
      ),
      scheme_code = peer,
      scheme_year = dplyr::if_else(
        stringr::str_detect(run_stage, "Final"),
        paste0(peer_year, "*"), peer_year
      ),
      # model run
      run_scenario = scenario,
      run_stage,
      # mitigators
      mitigator_code = `Mitigator code`,
      mitigator_name = `Strategy name`,
      mitigator_variable = strategy,
      mitigator_activity_type = `Activity type`,
      mitigator_type = `Mitigator type`,
      mitigator_subset = `Strategy subset`,
      mitigator_group = Grouping,
      # mitigator value selections
      value_lo = value_1,
      value_hi = value_2,
      value_mid = value_lo + (value_hi - value_lo) / 2,
      value_range = value_hi - value_lo,
      value_point_or_range = dplyr::if_else(
        value_hi - value_lo == 0,
        "point",
        "range"
      ),
      value_time_profile = time_profile,
      # years
      year_baseline = baseline_year,
      year_horizon = horizon_year,
      year_range = year_horizon - year_baseline,
      # national elicitation exercise
      nee_p10 = percentile10,
      nee_p90 = percentile90,
      nee_p50 = nee_p10 - (nee_p10 - nee_p90) / 2,
      nee_mean = mean
    )

  data_prepared |>
    dplyr::select(
      tidyselect::starts_with("scheme"),
      tidyselect::starts_with("run"),
      tidyselect::starts_with("mitigator"),
      tidyselect::starts_with("value"),
      tidyselect::starts_with("year"),
      tidyselect::starts_with("nee")
    ) |>
    dplyr::arrange(scheme_name, mitigator_code)

}

get_all_schemes <- function(dat) {
  dat |>
    shiny::req() |>
    dplyr::distinct(scheme_name, scheme_code) |>
    dplyr::filter(!is.na(scheme_code)) |>
    dplyr::mutate(scheme_name = glue::glue("{scheme_name} ({scheme_code})")) |>
    dplyr::arrange(scheme_name) |>
    tibble::deframe()
}

get_all_mitigators <- function(dat) {
  dat |>
    shiny::req() |>
    dplyr::distinct(mitigator_name, mitigator_code) |>
    dplyr::filter(!is.na(mitigator_code)) |>
    dplyr::mutate(
      mitigator_name = glue::glue("{mitigator_code}: {mitigator_name}")
    ) |>
    dplyr::arrange(mitigator_code) |>
    tibble::deframe()
}

get_all_activity_types <- function(dat) {
  dat |>
    shiny::req() |>
    dplyr::distinct(mitigator_activity_type) |>
    dplyr::pull() |>
    sort()
}

get_all_mitigator_types <- function(dat) {
  dat |>
    shiny::req() |>
    dplyr::distinct(mitigator_type) |>
    dplyr::pull() |>
    sort()
}

get_all_mitigator_subsets <- function(dat) {
  dat |>
    shiny::req() |>
    dplyr::distinct(mitigator_subset) |>
    dplyr::pull() |>
    sort()
}

#' Get a lookup table of participating Trusts
#'
#' Read a csv lookup file from Azure storage and do some clean-up to ensure
#' one row per Trust.
#'
#' @param container_support The Azure container for supporting information, as obtained by `get_container()` from `fct_tabulate.R`
#'
#' @return Tibble of data listing participating Trusts
#' @export
get_trust_lookup <- function(container_support) {

  trust_lookup <-
    # read the data from Azure
    AzureStor::storage_read_csv(
      container = container_support,
      file = "nhp-scheme-lookup.csv",
      show_col_types = FALSE
    ) |>
    # Imperial College (RYJ) appears three times due to different hospital
    # sites, so simplify to one row
    dplyr::mutate(
      `Name of Hospital site` = dplyr::case_match(
        `Trust ODS Code`,
        'RYJ' ~ 'Imperial',
        .default = `Name of Hospital site`
      )
    ) |>
    # Ensure one row per trust - deals with Hampshire which appears twice
    dplyr::distinct(`Trust ODS Code`, .keep_all = TRUE) |>
    # Sort
    dplyr::arrange(`Trust ODS Code`)

  return(trust_lookup)
}
