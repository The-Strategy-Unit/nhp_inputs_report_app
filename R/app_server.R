#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # Prepare data ----

  ## Read data ----

  container_results <-
    get_container(container_name = Sys.getenv("AZ_STORAGE_CONTAINER_RESULTS"))

  container_support <-
    get_container(container_name = Sys.getenv("AZ_STORAGE_CONTAINER_SUPPORT"))

  board <- pins::board_connect()
  params <- pins::pin_read(board, name = "matt.dray/nhp_tagged_runs_params")
  runs_meta <- pins::pin_read(board, name = "matt.dray/nhp_tagged_runs_meta")
  rates_data <- pins::pin_read(board, name = "thomas.jemmett/inputs_app_rates_data_v2-1")
  extracted_params <- extract_params(params, runs_meta)

  skeleton_table <- prepare_skeleton_table(extracted_params)

  trust_code_lookup <- get_trust_lookup(container_support = container_support)

  mitigator_lookup <- container_support |>
    AzureStor::storage_read_csv("mitigator-lookup.csv", show_col_types = FALSE)

  nee_results <- container_support |> read_nee("nee_table.rds")

  peers <- container_support |>
    AzureStor::storage_load_rds("trust-peers.rds") |>
    dplyr::rename(scheme = procode)

  # inputs app yml (contains meta data)
  yaml <- yaml::read_yaml(
    file = here::here('inst', 'ref_inputs-golem-config.yml'),
    eval.expr = FALSE
  )
  yaml_df <- get_mitigator_baseline_description(yaml)

  ## Prep data ----
  dat <- populate_table(
    skeleton_table,
    extracted_params,
    trust_code_lookup,
    mitigator_lookup,
    nee_results,
    yaml_df
  )

  all_schemes <- get_all_schemes(dat)


  # Reactive Values
  ra <- shiny::reactiveValues(
    heatmap_min_height = NA,
    pointrange_min_height = NA,
    mixturedist_min_height = NA
  )
  ra$heatmap_min_height <- 100
  ra$pointrange_min_height <- 100
  ra$mixturedist_min_height <- 100

  # Reactives ----

  # ensure dat reflects the user's preferred view
  dat_reactive <- shiny::reactive({
    dat_return <- update_dat_values(
      dat = dat,
      values_displayed = input$values_displayed,
      include_point_estimates = input$include_point_estimates
    )

    return(dat_return)
  })

  dat_mixture_distributions <- shiny::reactive({

    dat_return <- dat_reactive() |>
      get_mixture_distributions_dat()

    return(dat_return)
  })

  # pre-calculate mixture distributions (avoid re-calculating)
  #dat_mixture_distributions <- get_mixture_distributions_dat(dat = dat())

  peer_set <- shiny::reactive({
    peers |>
      dplyr::filter(scheme == input$focus_scheme) |>
      dplyr::pull(peer)
  })

  dat_filtered <- reactive({

    if (input$activity_type != "All") {
      dat_return <- dat_reactive() |>
        dplyr::filter(mitigator_activity_type == input$activity_type)
    } else {
      dat_return <- dat_reactive()
    }

    # standardise values to 2041 if requested
    if (input$standardise_2041) {
      dat_return <- dat_return |>
        dplyr::mutate(
          value_lo = forecast_value(
            year_baseline = year_baseline,
            year_horizon = year_horizon,
            year_forecast = 2041,
            value_horizon = value_lo,
            value_displayed = input$values_displayed
          ),
          value_hi = forecast_value(
            year_baseline = year_baseline,
            year_horizon = year_horizon,
            year_forecast = 2041,
            value_horizon = value_hi,
            value_displayed = input$values_displayed
          ),
          value_mid = (value_hi - ((value_hi - value_lo) / 2)) |>
            round(digits = 3)
        )
    }

    return(dat_return)

  })

  available_mitigators <- reactive({
    dat_filtered() |> get_all_mitigators()
  })

  available_mitigator_groups <- reactive({
    dat_filtered() |> get_all_mitigator_groups()
  })

  mitigator_group_set <- reactive({
    dat_filtered() |>
      dplyr::filter(mitigator_group == input$mitigator_groups) |>
      dplyr::distinct(mitigator_code) |>
      dplyr::pull()
  })

  ## dat_selected_pointrange ----
  dat_selected_pointrange <- shiny::reactive({

    # ensure at least one selection is made for scheme and mitigator
    shiny::validate(
      shiny::need(input$schemes, message = "Select at least one scheme.")
    )

    shiny::validate(
      shiny::need(input$mitigators, message = "Select at least one mitigator.")
    )

    # get filtered data
    dat <- dat_filtered()

    # plot scheme names in reverse order so displayed correctly in ggplot
    if (!input$toggle_invert_facets) {
      dat <- dat |> dplyr::mutate(scheme_name = forcats::fct_rev(scheme_name))
    }

    # plot mitigators in order
    if (input$toggle_invert_facets) {
      dat <- dat |>
        dplyr::mutate(
          mitigator_code = forcats::fct_rev(mitigator_code),
          mitigator_name = stats::reorder(mitigator_name, as.numeric(mitigator_code)) # order mitigator_name to match mitigator_code
        )
    }

    # further filtering
    dat <- dat |>
      dplyr::filter(
        scheme_code %in% input$schemes,
        mitigator_code %in% input$mitigators
      )

    # show an aggregate summary where requested
    if (input$toggle_aggregate_summary) {

      # summarise all-but-the-focal scheme
      if (input$toggle_aggregate_summary_minmaxrange) {

        # use the full min/max range
        dat_summary <-
          dat |>
          dplyr::filter(
            !scheme_code %in% input$focus_scheme, # exclude focal scheme
            !is.na(value_mid) # need at least a mid-point
          ) |>
          dplyr::summarise(
            scheme_code = 'All',
            scheme_name = 'Summary ●',
            value_mid = mean(value_mid, na.rm = TRUE),
            value_lo = min(value_lo, na.rm = TRUE),
            value_hi = max(value_hi, na.rm = TRUE),
            .by = c(mitigator_code, mitigator_name)
          )

      } else {

        # use average values
        dat_summary <-
          dat |>
          dplyr::filter(
            !scheme_code %in% input$focus_scheme, # exclude focal scheme
            !is.na(value_mid) # need at least a mid-point
          ) |>
          dplyr::summarise(
            scheme_code = 'All',
            scheme_name = 'Summary ●',
            value_mid = mean(value_mid, na.rm = TRUE),
            value_lo = mean(value_lo, na.rm = TRUE),
            value_hi = mean(value_hi, na.rm = TRUE),
            .by = c(mitigator_code, mitigator_name)
          )
      }

      # add in nee data
      dat_summary <-
        dat_summary |>
        dplyr::left_join(
          y = dat |>
            dplyr::select(mitigator_code, nee_mean, nee_p90, nee_p10) |>
            dplyr::distinct(),
          by = 'mitigator_code'
        )

      # add summary to dat as a new row
      dat <- dplyr::bind_rows(
        dat,
        dat_summary
      ) |>
        # sort schemes by name alphabetically and the summary placed last
        dplyr::mutate(
          scheme_name = scheme_name |>
            base::factor() |>
            forcats::fct_relevel('Summary ●', after = Inf) |>
            forcats::fct_rev(),

          # sort scheme codes to match scheme names
          scheme_code = scheme_code |>
            base::factor(levels = unique(scheme_code[order(scheme_name)]))
        )
    }

    # add colour coding
    dat <-
      dat |>
      dplyr::mutate(
        point_colour = dplyr::case_when(
          scheme_code == input$focus_scheme ~ 'red',
          scheme_code == 'All' ~ 'blue',
          .default = 'black'
        )
      )

    # return the result
    return(dat)

  })

  ## dat_selected_heatmap ----
  dat_selected_heatmap <- shiny::reactive({

    shiny::validate(
      shiny::need(input$schemes, message = "Select at least one scheme.")
    )

    shiny::validate(
      shiny::need(input$mitigators, message = "Select at least one mitigator.")
    )

    dat <- dat_filtered()
    dat <-
      dat |>
      prepare_heatmap_dat(
        dat = _,
        mitigator_codes = input$mitigators,
        scheme_codes = input$schemes,
        focal_scheme_code = input$focus_scheme,
        heatmap_type = input$heatmap_type,
        scheme_order = input$heatmap_scheme_order,
        mitigator_order = input$heatmap_mitigator_order,
        values_displayed = input$values_displayed,
        toggle_heatmap_nee = input$toggle_heatmap_nee,
        toggle_heatmap_aggregate_summaries = input$toggle_heatmap_aggregate_summaries
      )

    return(dat)

  })


  dat_selected_mixture_distributions <- shiny::reactive({

    shiny::validate(
      shiny::need(input$schemes, message = "Select at least one scheme.")
    )

    shiny::validate(
      shiny::need(input$mitigators, message = "Select at least one mitigator.")
    )

    # dat <- dat_filtered() |>
    #   dplyr::filter(
    #     mitigator_code %in% input$mitigators
    #   )
    #
    # dat <- get_mixture_distributions_dat(dat = dat)

    # using pre-calculated mixture distributions, filtered for selected mitigators
    dat <- dat_mixture_distributions() |>
      dplyr::filter(
        mitigator_code %in% input$mitigators
      )

  })


  # Observers ----

  ## Updates ----

  shiny::observe({
    shiny::updateSelectInput(
      session,
      "focus_scheme",
      choices = all_schemes,
      selected = all_schemes[1]
    )
  })

  shiny::observe({

    scheme_and_peers <- c(input$focus_scheme, peer_set())
    selected_schemes <- all_schemes[which(all_schemes %in% scheme_and_peers)] |>
      tibble::enframe() |>
      dplyr::arrange(name) |>  # sort on name, not code
      tibble::deframe()

    if (input$toggle_all_schemes) selected_schemes <- all_schemes

    shiny::updateSelectInput(
      session,
      "schemes",
      choices = all_schemes,
      selected = selected_schemes
    )

  })

  shiny::observe({
    shiny::updateSelectInput(
      session,
      "mitigator_groups",
      choices = available_mitigator_groups(),
      selected = available_mitigator_groups()[1]
    )
  })

  shiny::observe({
    shiny::updateSelectInput(
      session,
      "mitigators",
      choices = available_mitigators(),
      selected = mitigator_group_set()
    )
  })

  ## Enablers ----

  shiny::observe({

    if (input$heatmap_type == "value_binary") {
      # enable

      # disable
      shinyjs::disable("toggle_horizon_heatmap")
      shinyjs::disable("toggle_heatmap_scale_fill_by_mitigator")
      shinyjs::disable("toggle_heatmap_nee")
      shinyjs::disable("toggle_heatmap_aggregate_summaries")

      # show
      shinyjs::show("heatmap_binary_colour")

      # hide
      shinyjs::hide("heatmap_value_colour_low")
      shinyjs::hide("heatmap_value_colour_high")
    }

    if (input$heatmap_type != "value_binary") {
      # enable
      shinyjs::enable("toggle_horizon_heatmap")
      shinyjs::enable("toggle_heatmap_scale_fill_by_mitigator")
      shinyjs::enable("toggle_heatmap_nee")
      shinyjs::enable("toggle_heatmap_aggregate_summaries")

      # disable

      # show
      shinyjs::show("heatmap_value_colour_low")
      shinyjs::show("heatmap_value_colour_high")

      # hide
      shinyjs::hide("heatmap_binary_colour")
    }

    # disable 'summary full range' switch if 'summary' is disabled
    if (input$toggle_aggregate_summary) {
      shinyjs::enable("toggle_aggregate_summary_minmaxrange")
    } else {
      shinyjs::disable("toggle_aggregate_summary_minmaxrange")
    }

  })

  # Renders ----

  ## Plots ----

  ### pointrange ----

  # adjust the 'number of facet columns' slider to match the number of facets
  # i.e. so can select up to a minimum of one facet per row or list all on a
  # single row
  shiny::observe({

    dat <- dat_selected_pointrange()

    # what is the maximum number of columns to facet over?
    max_facet_cols <- if (input$toggle_invert_facets) {
      # count the number of selected schemes
      dat |>
        dplyr::filter(scheme_code %in% input$schemes) |>
        dplyr::pull(scheme_code) |>
        unique() |>
        length()

    } else {
      # count the number of selected mitigators
      dat |>
        dplyr::filter(mitigator_code %in% input$mitigators) |>
        dplyr::pull(mitigator_code) |>
        unique() |>
        length()
    }

    # update the ui
    shiny::updateSliderInput(
      inputId = 'facet_columns',
      max = max_facet_cols
    )

  })

  # adjust pointrange height in response to the number of facet rows
  # this will be base_height + (y_rows * facet_rows)
  shiny::observe({

    dat <- dat_selected_pointrange()

    # how many columns are specified?
    facet_cols <- input$facet_columns

    # how many mitigators are shown
    n_mitigators <- dat |>
      dplyr::filter(mitigator_code %in% input$mitigators) |>
      dplyr::pull(mitigator_code) |>
      unique() |>
      length()

    # how many schemes are shown
    n_schemes <- dat |>
      dplyr::filter(scheme_code %in% input$schemes) |>
      dplyr::pull(scheme_code) |>
      unique() |>
      length()

    # determine the number of rows and facet rows
    if (input$toggle_invert_facets) {
      # facetted on scheme with mitigator in the y-axis

      n_rows <- n_mitigators
      n_facets <- n_schemes

    } else {
      # facetted on mitigator with schemes in the y-axis

      n_rows <- n_schemes
      n_facets <- n_mitigators

    }

    # estimate facet rows
    n_facet_rows <- ceiling(n_facets / facet_cols)

    # specify the height for the pointrange plot
    base_height <- 150

    # update reactive values
    # update the min_height to base + pro-rata for each row and facet row
    ra$pointrange_min_height <- base_height +
      (n_facet_rows * n_rows * 35) + # add room for each row on each facet
      (n_facet_rows * 100) # add room for facet labels on each row

  })

  # wrap the plot call in an observer to enable the dynamic height setting
  shiny::observe({
    output$pointrange <- shiny::renderPlot({

      shiny::validate(
        shiny::need(
          nrow(dat_selected_pointrange()) > 0,
          message = "Insufficient data for this plot."
        )
      )

      dat_selected_pointrange() |> plot_pointrange(input)

    }, height = ra$pointrange_min_height)
  })


  ### heatmap ----

  # adjust heatmap height in response to the number of mitigators to display
  shiny::observe({

    dat <- dat_selected_heatmap()

    # count the number of selected mitigators
    temp_mitigator_count <- dat |>
      dplyr::filter(mitigator_code %in% input$mitigators) |>
      dplyr::pull(mitigator_code) |>
      unique() |>
      length()

    # count the number of selected schemes
    temp_scheme_count <- dat |>
      dplyr::filter(scheme_code %in% input$schemes) |>
      dplyr::pull(scheme_code) |>
      unique() |>
      length()

    # decide the base height (the space taken up by scheme names and legend)
    # NB, scheme names are rotated as the scheme count increases to fit them
    # on screen, so they take up more vertical space
    if (temp_scheme_count < 6) {
      base_height <- 100
    } else if (temp_scheme_count < 12) {
      base_height <- 200
    } else {
      base_height <- 300
    }

    # update reactive values
    # update the min_height to base + pro-rata for each mitigator
    ra$heatmap_min_height <- base_height + (temp_mitigator_count * 55)
  })

  # produce the heatmap plot
  output$heatmap <- plotly::renderPlotly({

    shiny::validate(
      shiny::need(
        nrow(dat_selected_heatmap()) > 0,
        message = "Insufficient data for this plot."
      )
    )

    dat_selected_heatmap() |>
      plot_heatmap(
        # data
        focal_scheme_code = input$focus_scheme,

        # options
        toggle_mitigator_name = input$toggle_mitigator_name,
        toggle_scale_fill_by_mitigator = input$toggle_heatmap_scale_fill_by_mitigator,
        values_displayed = input$values_displayed,
        heatmap_type = input$heatmap_type,

        # formatting
        colour_binary = input$heatmap_binary_colour,
        colour_value_low = input$heatmap_value_colour_low,
        colour_value_high = input$heatmap_value_colour_high,
        plot_height = ra$heatmap_min_height
      )

  })

  ### density functions ----
  # adjust distribution height in response to the number of mitigators to display
  shiny::observe({

    #dat <- dat_selected_mixture_distributions()
    dat <- dat_selected_heatmap()

    # count the number of selected mitigators
    temp_mitigator_count <- dat |>
      dplyr::filter(mitigator_code %in% input$mitigators) |>
      dplyr::pull(mitigator_code) |>
      unique() |>
      length()

    # update reactive values
    # update the min_height to base + pro-rata for each mitigator
    ra$mixturedist_min_height <- (temp_mitigator_count * 200)
  })

  # wrap the plot call in an observer to enable the dynamic height setting
  shiny::observe({
    output$mixture_distributions <- shiny::renderPlot({
      plot_mixture_distributions(
        dat_selected_mixture_distributions = dat_selected_mixture_distributions(),
        dat_filtered = dat_filtered(),
        dat_focal_scheme_code = input$focus_scheme,
        input = input
      )
    }, height = ra$mixturedist_min_height)
  })

  # output$mixture_distributions <- shiny::renderPlot({
  #   dat_selected_mixture_distributions() |>
  #     plot_mixture_distributions(input)
  # }, height = 10000)

  ### contextual baseline -----

  output$contextual_baseline <- plotly::renderPlotly({

    # notify user if no data is available
    shiny::validate(
      shiny::need(input$schemes, message = "Select at least one scheme.")
    )

    shiny::validate(
      shiny::need(input$mitigators, message = "Select at least one mitigator.")
    )

    # plot the baseline contextual
    dat_filtered() |>
      plot_baseline_comparison(
        rates_data = rates_data,
        mitigator_codes = input$mitigators,
        focal_scheme_code = input$focus_scheme,
        rate_title = "Baseline rate",
        value_title = input$values_displayed,
        trendline = FALSE, # trendlines don't seem to work
        range = input$toggle_contextual_baseline_range,
        scheme_label = input$toggle_contextual_baseline_schemecode,
        quadrants = input$toggle_contextual_baseline_quadrants,
        facet_columns = 1, # not enabling this feature - feel it is best for one row per mitigator
        facet_height_px = input$slider_contextual_baseline_height
      )

  })

  ### contextual trendline ----
  output$contextual_trendline <- plotly::renderPlotly({

    # notify user if no data is available
    shiny::validate(
      shiny::need(input$schemes, message = "Select at least one scheme.")
    )

    shiny::validate(
      shiny::need(input$mitigators, message = "Select at least one mitigator.")
    )

    # plot the trendline contextual
    dat_filtered() |>
      dplyr::filter(mitigator_code %in% input$mitigators) |>
      plot_faceted_trendlines(
        rates_data = rates_data,
        mitigator_codes = input$mitigators,
        focal_scheme_code = input$focus_scheme,
        scheme_codes = input$schemes,
        show_other_schemes = input$toggle_contextual_trendline_otherschemes,
        show_horizon_timeline = input$toggle_contextual_trendline_horizon_timeline,
        show_horizon_overlay = input$toggle_contextual_trendline_horizon_overlay,
        show_prebaseline_average = input$toggle_contextual_trendline_average,
        facet_height_px = input$slider_contextual_trendline_height
      )

  })

  ## Tables ----

  output$raw_data_dt <- DT::renderDT({
    # NB, using the original dat tibble as want to present the un-edited data
    dat |> make_raw_dt()
  }, server = FALSE) # to download all data via CSV button

  output$mitigator_lookup_dt <- DT::renderDT({
    mitigator_lookup |> make_mitigator_dt()
  }, server = FALSE) # to download all data via CSV button

  output$scheme_lookup_dt <- DT::renderDT({
    trust_code_lookup |> make_scheme_dt()
  }, server = FALSE) # to download all data via CSV button

  output$mitigator_uptake_dt <- DT::renderDT({
    make_mitigator_uptake_dt(
      dat = dat,
      selected_schemes = input$schemes
    )
  })

  output$scheme_uptake_dt <- DT::renderDT({
    make_scheme_uptake_dt(
      dat = dat,
      selected_mitigators = input$mitigators,
      selected_schemes = input$schemes,
      focal_scheme = input$focus_scheme
    )
  })

}
