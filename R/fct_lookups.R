make_raw_dt <- function(dat) {

  dat_prepared <- dat |>
    dplyr::filter(!is.na(value_lo)) |>  # only want mitigators selected by schemes
    dplyr::mutate(
      # remove any pencil emojis from scheme names
      scheme_name = stringr::str_remove(
        string = scheme_name,
        pattern = ' ✏️'
      ),
      dplyr::across(
        c(
          tidyselect::starts_with("scheme"),
          tidyselect::starts_with("run"),
          tidyselect::starts_with("mitigator"),
          value_point_or_range,
          value_time_profile
        ),
        factor  # enables discrete selection in interactive table
      ),
    )

  dat_prepared |>
    DT::datatable(
      fillContainer = TRUE,
      filter = "top",
      rownames = FALSE,
      selection = "none",
      extensions = "Buttons",
      options = list(
        dom = "Bftp",
        buttons = list(
          list(
            extend = "csv",
            filename = paste0(Sys.Date(), "_mitigator-comparison-data"),
            text = "Download (CSV)"
          )
        )
      ),
      callback = DT::JS(
        "$('button.buttons-csv').css('background','#337ab7');
        $('button.buttons-csv').css('border','#2e6da4');
        return table;"
      )
    ) |>
    DT::formatRound(  # present to 3 dp, download will have unrestricted dp
      columns = c(
        "value_lo", "value_hi", "value_mid", "value_range",
        "nee_p10", "nee_p90", "nee_p50", "nee_mean"
      ),
      digits = 3
    )

}

make_mitigator_dt <- function(mitigator_lookup) {

  mitigator_lookup_prepared <- mitigator_lookup |>
    dplyr::select(
      `Mitigator code`,
      `Mitigator name` = `Strategy name`,
      `Mitigator variable` = `Strategy variable`,
      `Mitigator type`,
      `Activity type`,
      `Mitigator subset` = `Strategy subset`,
      `Mitigator grouping` = Grouping
    ) |>
    dplyr::mutate(dplyr::across(tidyselect::everything(), factor))

  mitigator_lookup_prepared |>
    DT::datatable(
      fillContainer = TRUE,
      filter = "top",
      rownames = FALSE,
      selection = "none",
      extensions = "Buttons",
      options = list(
        dom = "Bftp",
        buttons = list(
          list(
            extend = "csv",
            filename = paste0(Sys.Date(), "_mitigator-lookup"),
            text = "Download (CSV)"
          )
        )
      ),
      callback = DT::JS(
        "$('button.buttons-csv').css('background','#337ab7');
        $('button.buttons-csv').css('border','#2e6da4');
        return table;"
      )
    )

}


#' Make the scheme lookup DT object
#'
#' Renders a DT object listing schemes represented in the app with filters and
#' a CSV download button.
#'
#' @param trust_code_lookup Tibble of scheme data - as provided by the `get_trust_lookup()` function in `fct_tabulate.R`
#'
#' @return DT object listing schemes
#' @export
make_scheme_dt <- function(trust_code_lookup) {

  schemes_prepared <- trust_code_lookup |>
    dplyr::mutate(
      `Scheme name` = `Name of Hospital site`,
      `Trust name` = `Name of Trust`,
      `Scheme code` = `Trust ODS Code`,
      .keep = "none"
    ) |>
    dplyr::mutate(across(tidyselect::everything(), factor)) |>
    dplyr::arrange(`Scheme name`)

  schemes_prepared |>
    DT::datatable(
      fillContainer = TRUE,
      filter = "top",
      rownames = FALSE,
      selection = "none",
      extensions = "Buttons",
      options = list(
        pageLength = 100, # show all schemes
        dom = "Bft", # excluding the paging buttons
        buttons = list(
          list(
            extend = "csv",
            filename = paste0(Sys.Date(), "_mitigator-comparison-data"),
            text = "Download (CSV)"
          )
        )
      ),
      callback = DT::JS(
        "$('button.buttons-csv').css('background','#337ab7');
        $('button.buttons-csv').css('border','#2e6da4');
        return table;"
      )
    )

}


#' Make the scheme uptake DT object
#'
#' Renders a DT object showing the proportion of mitigators in use by each scheme.
#' Two rates are shown:
#' 1. covers all available mitigators,
#' 2. covers the subset of mitigators selected by the user.
#'
#' @param dat Tibble - the full prepared dataset for this app
#' @param selected_schemes Character vector - a list of mitigator_codes selected by the user
#'
#' @return DT object listing schemes and the proportions of mitigators in use by them
#' @export
make_mitigator_uptake_dt <- function(dat, selected_schemes) {

  dat |>
    # remove ampersand from mitigator names - causes issues with DT filters
    dplyr::mutate(
      mitigator_name = gsub(
        pattern = 'A&E',
        x = mitigator_name,
        replacement = 'ED'
      )
    ) |>
    # count schemes per mitigator
    dplyr::summarise(
      n_schemes_using_all = dplyr::n_distinct(
        scheme_code,
        na.rm = TRUE
      ),
      n_schemes_using_selected = dplyr::n_distinct(
        scheme_code[scheme_code %in% selected_schemes],
        na.rm = TRUE
      ),
      .by = c(mitigator_activity_type, mitigator_subset, mitigator_name)
    ) |>
    # convert to rate
    dplyr::mutate(
      # get denominators
      n_schemes_all = dplyr::n_distinct(
        dat$scheme_code,
        na.rm = TRUE
      ),
      n_schemes_selected = dplyr::n_distinct(
        dat$scheme_code[dat$scheme_code %in% selected_schemes],
        na.rm = TRUE
      ),

      # convert to rates
      n_schemes_using_all_rate = n_schemes_using_all / n_schemes_all,
      n_schemes_using_selected_rate = n_schemes_using_selected / n_schemes_selected
    ) |>
    # prepare for display
    dplyr::select(
      -c(n_schemes_using_all, n_schemes_all,
         n_schemes_using_selected, n_schemes_selected)
    ) |>
    dplyr::mutate(
      # convert mitigators to factors for drop-down selectors in DT
      mitigator_activity_type = mitigator_activity_type |> factor(),
      mitigator_subset = mitigator_subset |> factor(),
      mitigator_name = mitigator_name |> factor()
    ) |>
    # display as DT
    DT::datatable(
      rownames = FALSE,
      options = list(pageLength = 100, dom = 'ft'),
      fillContainer = TRUE,
      escape = TRUE,
      filter = 'top',
      colnames = c(
        'Activity type', 'Mitigator group', 'Mitigator',
        'Coverage (all schemes)', 'Coverage (selected schemes)'
      )
    ) |>
    DT::formatPercentage(
      columns = c('n_schemes_using_all_rate', 'n_schemes_using_selected_rate')
    )
}


#' Make the mitigator uptake DT object
#'
#' Renders a DT object showing the proportion of mitigators used by each scheme.
#' Two rates are shown:
#' 1. covers all available mitigators,
#' 2. covers the subset of mitigators selected by the user.
#'
#' @param dat Tibble - the full prepared dataset for this app
#' @param selected_schemes Character vector - a list of scheme_codes selected by the user
#' @param focal_scheme Character vector - the focal scheme_code
#'
#' @return DT object listing schemes and the proportions of mitigators in use by them
#' @export
make_scheme_uptake_dt <- function(dat, selected_mitigators, selected_schemes, focal_scheme) {

  dat |>
    # count schemes per mitigator
    dplyr::summarise(
      n_mitigators_using_all = dplyr::n_distinct(
        mitigator_code,
        na.rm = T
      ),
      n_mitigators_using_selected = dplyr::n_distinct(
        mitigator_code[mitigator_code %in% selected_mitigators],
        na.rm = T
      ),
      .by = c(scheme_code, scheme_name)
    ) |>
    # convert to rate
    dplyr::mutate(
      # get denominators
      n_mitigators_all = dplyr::n_distinct(
        dat$mitigator_code,
        na.rm = T
      ),
      n_mitigators_selected = dplyr::n_distinct(
        dat$mitigator_code[dat$mitigator_code %in% selected_mitigators],
        na.rm = T
      ),

      # convert to rates
      n_mitigators_using_all_rate = n_mitigators_using_all / n_mitigators_all,
      n_mitigators_using_selected_rate = n_mitigators_using_selected / n_mitigators_selected
    ) |>
    # prepare for display
    dplyr::select(
      -c(n_mitigators_using_all, n_mitigators_all,
         n_mitigators_using_selected, n_mitigators_selected)
    ) |>
    dplyr::mutate(
      # convert scheme details to factors for drop-down selectors in DT
      scheme_code = scheme_code |> factor(),
      scheme_name = scheme_name |> factor()
    ) |>
    dplyr::filter(!is.na(scheme_code)) |>
    # display as DT
    DT::datatable(
      rownames = FALSE,
      options = list(pageLength = 100, dom = 'Bft'),
      fillContainer = TRUE,
      escape = TRUE,
      style = 'default', # needed to ensure formatStyle works as expected - due to clashes with bslib & bootstrap theme
      colnames = c(
        'Scheme code', 'Scheme name', 'Coverage (all mitigators)',
        'Coverage (selected mitigators)'
      ),
      filter = 'top'
    ) |>
    DT::formatPercentage(
      columns = c('n_mitigators_using_all_rate', 'n_mitigators_using_selected_rate')
    ) |>
    # style selected schemes in bold
    DT::formatStyle(
      columns = 'scheme_code',
      target = 'row',
      # highlight all selected schemes in bold
      fontWeight = DT::styleEqual(levels = c(selected_schemes), 'bold', 'normal'),
      # highlight focal scheme in red too
      color = DT::styleEqual(levels = focal_scheme, 'red', 'black')
    )
}
