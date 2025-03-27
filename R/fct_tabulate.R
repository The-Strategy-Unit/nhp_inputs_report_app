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
#' @param yaml_df Tibble - output of `get_mitigator_baseline_description()`
#'
#' @return Tibble of data
#' @export
populate_table <- function(
    skeleton_table,
    extracted_params,
    trust_code_lookup,
    mitigator_lookup,
    nee_results,
    yaml_df
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
    ) |>
    dplyr::left_join(
      y = yaml_df |>
        dplyr::select(strategy_variable, mitigator_activity_title = y_axis_title) |>
        dplyr::distinct(),
      by = dplyr::join_by(strategy == strategy_variable)
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
      mitigator_group = Grouping,
      mitigator_activity_title = mitigator_activity_title,
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

  # prepare the percent mitigated (pm) and prediction interval (pi) redundant columns
  data_prepared <-
    data_prepared |>
    dplyr::mutate(
      # prediction_interval
      pi_value_lo = value_lo,
      pi_value_hi = value_hi,
      pi_value_mid = value_mid,
      pi_nee_p10 = round(nee_p10, 3),
      pi_nee_p90 = round(nee_p90, 3),
      pi_nee_p50 = round(nee_p50, 3),
      pi_nee_mean = round(nee_mean, 3),

      # percent mitigated
      pm_value_lo = round(1 - value_lo, 3),
      pm_value_hi = round(1 - value_hi, 3),
      pm_value_mid = round(1 - value_mid, 3),
      pm_nee_p10 = round(1 - nee_p10, 3),
      pm_nee_p90 = round(1 - nee_p90, 3),
      pm_nee_p50 = round(1 - nee_p50, 3),
      pm_nee_mean = round(1 - nee_mean, 3)
    )

  # organise ready for output
  data_prepared |>
    dplyr::select(
      tidyselect::starts_with("scheme"),
      tidyselect::starts_with("run"),
      tidyselect::starts_with("mitigator"),
      tidyselect::starts_with("value"),
      tidyselect::starts_with("year"),
      tidyselect::starts_with("nee"),
      tidyselect::starts_with("pi_"),
      tidyselect::starts_with("pm_")
    ) |>
    dplyr::arrange(scheme_name, mitigator_code)

}

#' Update dat to reflect the user's preferred values
#'
#' Users are able to select from 'Percent of activity mitigated' and '80%
#' prediction interval' in the global settings sidebar.
#'
#' This function updates dat to reflect the user's preferred values.
#'
#' @param dat Tibble of data - as produced from `populate_table` in `fct_tabulate.R`
#' @param values_displayed Character - the value in `input$values_displayed` indicating the user's preferred view
#' @param include_point_estimates Boolean - the value in `input$include_point_estimates` indicating whether point-estimates of 0% mitigation (100% prediction) are to be included
#' @param focal_scheme_code Character - the scheme code of the focal scheme
#'
#' @return Tibble of data with the 'value_' and 'nee_' fields updated to match the user's preferred values
update_dat_values <- function(dat,
                              values_displayed,
                              include_point_estimates = FALSE,
                              focal_scheme_code) {
  if (values_displayed == "Percent of activity mitigated") {
    dat <-
      dat |>
      dplyr::mutate(
        value_lo = pm_value_hi, # low becomes high if viewing percent mitigated
        value_hi = pm_value_lo, # high becomes low if viewing percent mitigated
        value_mid = pm_value_mid,
        nee_p10 = pm_nee_p90, # low becomes high if viewing percent mitigated
        nee_p90 = pm_nee_p10, # high becomes low if viewing percent mitigated
        nee_p50 = pm_nee_p50,
        nee_mean = pm_nee_mean
      )
  } else {
    dat <-
      dat |>
      dplyr::mutate(
        value_lo = pi_value_lo,
        value_hi = pi_value_hi,
        value_mid = pi_value_mid,
        nee_p10 = pi_nee_p10,
        nee_p90 = pi_nee_p90,
        nee_p50 = pi_nee_p50,
        nee_mean = pi_nee_mean
      )
  }

  if (!include_point_estimates) {
    # we need to exclude point estimates of value = 1
    # NB, I don't know of a one-step way of filtering data to exclude records
    # where two or more conditions are met.
    # So, here we will do this in a two-step process, 1) define the records to
    # exclude then 2) filter out these records

    # Create a key so we can refer to specific rows
    dat_temp <-
      dat |>
      dplyr::mutate(key = dplyr::row_number())

    # Define which records meet the definition of a point-estimate
    dat_exclude <-
      dat_temp |>
      dplyr::filter(pi_value_lo == 1 & pi_value_hi == 1)

    # Define `Dat` to exclude point-estimates
    dat <-
      dat_temp |>
      dplyr::filter(!key %in% dat_exclude$key) |>
      dplyr::select(-key)
  }

  # identify the focal scheme
  dat <-
    dat |>
    dplyr::mutate(
      scheme_name = dplyr::if_else(
        condition = scheme_code %in% dplyr::coalesce(focal_scheme_code, ""),
        true = glue::glue("{scheme_name} 📌"),
        false = scheme_name
      )
    )

  return(dat)
}

#' Forecast a scheme's value for a given forecast year
#'
#' Uses a simple linear model of a scheme's value - either `value_hi` or
#' `value_lo` to predict the value at `year_forecast`.
#'
#' Many schemes have a horizon year of 2041, but some have other horizons, such
#' as 2030 or 2035 which can make it tricky to compare schemes.
#'
#' This function predicts what the value would be at 2041, providing a
#' standardised value that makes comparing schemes a little easier.
#'
#' @param year_baseline The year the activity starts from, defaults to 2019
#' @param year_horizon The year the value is predicted for, often 2041 but can be different
#' @param year_forecast The year to forecast for, defaults to 2041
#' @param value_horizon The value at the horizon year
#' @param value_displayed The content of `input$values_displayed` - whether values are displayed as prediction intervals or percent mitigated, as these require differnet calculation approaches
#'
#' @return numeric predicted value between 0 and 1 rounded to 3 decimal places at `year_forecast`
forecast_value <- function(
    year_baseline = 2019,
    year_horizon,
    year_forecast = 2041,
    value_horizon,
    value_displayed
) {

  # gather some details
  horizon_length <- (year_horizon - year_baseline)
  forecast_length <- (year_forecast - year_baseline)

  # calculate linear prediction for year_forecast - adjust to whether pm or pi in use
  if (value_displayed == 'Percent of activity mitigated') {
    value_forecast <- ((value_horizon / horizon_length) * forecast_length)
  } else {
    value_forecast <- 1 - (((1 - value_horizon) / horizon_length) * forecast_length)
  }

  # bound value to within 0 or 1
  value_forecast <- dplyr::case_when(
    value_forecast > 1 ~ 1,
    value_forecast < 0 ~ 0,
    .default = value_forecast
  )

  # limit to 3 d.p.
  value_forecast <- round(value_forecast, digits = 3)

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

get_all_mitigator_groups <- function(dat) {
  dat |>
    shiny::req() |>
    dplyr::distinct(mitigator_group) |>
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

#' Get mitigator baseline descriptions
#'
#' The baseline activity for each mitigator is not specified in the `dat`, but
#' is implicit. This function extracts the activity description, e.g.
#' 'Admissions per 1,000 population' or '% of Appointments that are Face-to-Face'
#' from the inputs app yaml file.
#'
#' Using Gabriel's code from here:#'
#' https://github.com/The-Strategy-Unit/nhp_schemes_report/blob/55de2b1e67394a74eeef1e63ae3b88720c281c3d/R/mitigator%20credibility%20from%20historical%20data.R#L101-L137
#'
#' @param yaml List - the yaml file used in the inputs app
#'
#' @returns Tibble of mitigator details including the `y_axis_title` activity
get_mitigator_baseline_description <- function(yaml) {

  # abbreviate the path to the mitigators_config list
  mitigator_yaml <- yaml$default$mitigators_config

  # map over the categories and create a data frame
  df_return <- purrr::map_dfr(
    .x = names(mitigator_yaml),
    function(category_name) {

      # add to a data frame with category, element, and y_axis title
      return(
        tibble::tibble(
          activity_type = mitigator_yaml[[category_name]]$activity_type,
          mitigator_type = mitigator_yaml[[category_name]]$mitigators_type,
          category = category_name,
          strategy_variable = names(mitigator_yaml[[category_name]]$strategy_subset),
          y_axis_title = mitigator_yaml[[category_name]]$y_axis_title
        )
      )
    }
  )

  # return the result
  return(df_return)
}


#' Prepare Mitigator lookup
#'
#' Prepare the mitigator lookup file for subsequent use by
#' a) cleaning variable names
#' b) ordering the variables by cardinality
#'
#' @param mitigator_lookup Tibble of mitigator lookup data
#'
#' @return Tibble of mitigator lookups
prepare_mitigators <- function(mitigator_lookup) {

  # clean up the variable names
  mitigator_lookup <-
    mitigator_lookup |>
    janitor::clean_names()

  # find out the cardinality of each variable
  mit_order <-
    mitigator_lookup |>
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        \(.x) length(unique(.x))
      )
    ) |>
    tidyr::pivot_longer(cols = dplyr::everything()) |>
    dplyr::arrange(value)

  # order the variables by cardinality (less to more)
  mitigator_lookup <-
    mitigator_lookup |>
    dplyr::select(dplyr::any_of(mit_order$name)) |>
    dplyr::rename(mitigator_name = strategy_name) |>
    dplyr::mutate(mitigator_name = glue::glue('{mitigator_name} [{mitigator_code}]')) |>
    dplyr::select(-strategy_variable)

  return(mitigator_lookup)
}

#' Add to the selected mitigator list
#'
#' Adds newly selected mitigators to the list of currently selected mitigators
#'
#' @param df Tibble of mitigator reference data - e.g. `mitigator_reference`
#' @param selected_currently A list of mitigator_codes which are currently selected - as provided by `input$mitigators`
#' @param new_selections A list of newly selected mitigatoe
#'
#' @return Vector - mitigator labels with mitigator codes as the value
add_to_selected_mitigators <- function(df, selected_currently, new_selections) {

  df |>
    dplyr::filter(mitigator_code %in% c(selected_currently, new_selections)) |>
    dplyr::select(mitigator_name, mitigator_code) |>
    tibble::deframe()

}
