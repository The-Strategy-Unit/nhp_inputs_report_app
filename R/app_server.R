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
  extracted_params <- extract_params(params, runs_meta)

  skeleton_table <- prepare_skeleton_table(extracted_params)

  trust_code_lookup <- get_trust_lookup(container_support = container_support)

  mitigator_lookup <- container_support |>
    AzureStor::storage_read_csv("mitigator-lookup.csv", show_col_types = FALSE)

  nee_results <- container_support |> read_nee("nee_table.rds")

  peers <- container_support |>
    AzureStor::storage_load_rds("trust-peers.rds") |>
    dplyr::rename(scheme = procode)

  ## Prep data ----

  dat <- populate_table(
    skeleton_table,
    extracted_params,
    trust_code_lookup,
    mitigator_lookup,
    nee_results
  )

  all_schemes <- get_all_schemes(dat)

  # Reactive Values
  ra <- shiny::reactiveValues(
    heatmap_min_height = NA,
    pointrange_min_height = NA
  )
  ra$heatmap_min_height <- 100
  ra$pointrange_min_height <- 100

  # Reactives ----

  # The peers associated with the focal scheme
  peer_set <- shiny::reactive({
    peers |>
      dplyr::filter(scheme == input$focus_scheme) |>
      dplyr::pull(peer)
  })

  # The full set of available activity types in the raw data
  available_activity_types <- reactive(get_all_activity_types(dat))

  # The set of available mitigator types given the activity type
  mitigator_type_set <- reactive({
    dat |>
      dplyr::filter(mitigator_activity_type %in% input$activity_type) |>
      dplyr::distinct(mitigator_type) |>
      dplyr::pull()
  })

  # The set of available mitigator subsets given the activity and mitigator type
  mitigator_subsets_set <- reactive({
    dat |>
      dplyr::filter(
        mitigator_activity_type %in% input$activity_type,
        mitigator_type %in% input$mitigator_type,
      ) |>
      dplyr::distinct(mitigator_subset) |>
      dplyr::pull()
  })

  # The set of available mitigators given the activity and mitigator types and
  # mitigator subsets
  mitigator_set <- reactive({
    dat |>
      dplyr::filter(
        mitigator_activity_type %in% input$activity_type,
        mitigator_type %in% input$mitigator_type,
        mitigator_subset %in% input$mitigator_subset
      ) |>
      get_all_mitigators()  # prepare named vector
  })

  # Current dataset given filters
  dat_filtered <- reactive({
    dat |> dplyr::filter(mitigator_code %in% mitigator_set())
  })

  # Plots

  dat_selected_pointrange <- shiny::reactive({

    shiny::validate(
      need(input$schemes, message = "Select at least one scheme.")
    )

    shiny::validate(
      need(input$mitigators, message = "Select at least one mitigator.")
    )

    dat <- dat_filtered()

    if (input$toggle_horizon_pointrange) {
      dat <- dat |>
        dplyr::mutate(
          dplyr::across(
            c(value_lo, value_hi, value_mid),
            \(x) x / year_range
          )
        )
    }

    if (!input$toggle_invert_facets) {
      dat <- dat |> dplyr::mutate(scheme_name = forcats::fct_rev(scheme_name))
    }

    if (input$toggle_invert_facets) {
      dat <- dat |>
        dplyr::mutate(
          mitigator_code = forcats::fct_rev(mitigator_code),
          mitigator_name = stats::reorder(mitigator_name, as.numeric(mitigator_code)) # order mitigator_name to match mitigator_code
        )
    }

    dat |>
      dplyr::filter(
        scheme_code %in% input$schemes,
        mitigator_code %in% input$mitigators
      ) |>
      dplyr::mutate(
        point_colour = dplyr::if_else(
          scheme_code == input$focus_scheme,
          TRUE,
          FALSE
        )
      )

  })

  dat_selected_heatmap <- shiny::reactive({

    shiny::validate(
      need(input$schemes, message = "Select at least one scheme.")
    )

    shiny::validate(
      need(input$mitigators, message = "Select at least one mitigator.")
    )

    dat <- dat_filtered()

    if (input$toggle_horizon_heatmap) {
      dat <- dat |>
        dplyr::mutate(
          dplyr::across(
            c(value_lo, value_hi, value_mid),
            \(x) x / year_range
          )
        )
    }

    dat |>
      dplyr::mutate(
        value_binary = dplyr::if_else(!is.na(value_lo), 1, 0),
        dplyr::across(
          c(value_lo, value_hi, value_mid, value_range),
          \(x) janitor::round_half_up(x, 3)
        ),
        mitigator_code = forcats::fct_rev(mitigator_code),
        mitigator_name = stats::reorder(mitigator_name, as.numeric(mitigator_code)) # order mitigator_name to match mitigator_code
      ) |>
      tidyr::pivot_longer(
        c(value_lo, value_hi, value_mid, value_range, value_binary),
        names_to = "value_type",
        values_to = "value"
      ) |>
      dplyr::filter(
        value_type == input$heatmap_type,
        scheme_code %in%  input$schemes,
        mitigator_code %in% input$mitigators,
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

    selected_schemes <- c(input$focus_scheme, peer_set()) |> sort()

    if (input$toggle_all_schemes) {
      selected_schemes <- all_schemes
    }

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
      "activity_type",
      choices = available_activity_types(),
      selected = available_activity_types()[1]
    )
  })

  shiny::observe({
    shiny::updateSelectInput(
      session,
      "mitigator_type",
      choices = mitigator_type_set(),
      selected = mitigator_type_set()[1]
    )
  })

  shiny::observe({
    shiny::updateSelectInput(
      session,
      "mitigator_subset",
      choices = mitigator_subsets_set(),
      selected = mitigator_subsets_set()[1]
    )
  })

  shiny::observe({
    shiny::updateSelectInput(
      session,
      "mitigators",
      choices = mitigator_set(),
      selected = mitigator_set()
    )
  })

  ## Enablers ----

  shiny::observe({

    if (input$heatmap_type == "value_binary") {
      shinyjs::disable("toggle_horizon_heatmap")
    }

    if (input$heatmap_type != "value_binary") {
      shinyjs::enable("toggle_horizon_heatmap")
    }

    # disable NEE reference range checkbox when viewing standardised horizon
    if (input$toggle_horizon_pointrange) {
      shinyjs::disable("toggle_nee_reference_range")
    } else {
      shinyjs::enable("toggle_nee_reference_range")
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
      base_height <- 160
    } else {
      base_height <- 200
    }

    # update reactive values
    # update the min_height to base + pro-rata for each mitigator
    ra$heatmap_min_height <- base_height + (temp_mitigator_count * 55)
  })

  # wrap the plot call in an observer to enable the dynamic height setting
  shiny::observe({
    output$heatmap <- shiny::renderPlot({
      dat_selected_heatmap() |> plot_heatmap(input)
    }, height = ra$heatmap_min_height)
  })

  ## Tables ----

  output$raw_data_dt <- DT::renderDT({
    dat_filtered() |> make_raw_dt()
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
