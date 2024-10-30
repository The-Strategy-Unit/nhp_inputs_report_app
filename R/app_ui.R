#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  shiny::tagList(
    # External resources
    golem_add_external_resources(),
    shinyjs::useShinyjs(),
    # Application UI logic
    bslib::page_navbar(
      id = "page_navbar",
      title = bslib::tooltip(
        trigger = list(
          "NHP mitigator comparisons",
          bsicons::bs_icon("exclamation-triangle")
        ),
        md_file_to_html("app", "text", "warning.md")
      ),
      header = bslib::card(
        fill = FALSE,
        bslib::card_header(
          class = "bg-warning",
          bsicons::bs_icon("exclamation-triangle"),
          "Warning"
        ),
        md_file_to_html("app", "text", "warning.md")
      ),
      ## sidebar ----
      sidebar = bslib::sidebar(
        id = "sidebar",
        width = 400,
        bslib::accordion(
          id = "global_accordion",
          open = FALSE,
          ### scheme select ----
          bslib::accordion_panel(
            id = "accordion_schemes",
            title = "Select schemes",
            icon = bsicons::bs_icon("hospital"),
            shiny::selectInput(
              inputId = "focus_scheme",
              label = bslib::tooltip(
                trigger = list(
                  "Focal scheme",
                  bsicons::bs_icon("info-circle")
                ),
                "The scheme to highlight in plots."
              ),
              choices = NULL,
              selected = NULL,
              multiple = FALSE
            ),
            shiny::selectInput(
              inputId = "schemes",
              label = bslib::tooltip(
                trigger = list(
                  "Schemes to visualise",
                  bsicons::bs_icon("info-circle")
                ),
                "Defaults to peers of the selected scheme if 'Select all schemes?' is unchecked."
              ),
              choices = NULL,
              selected = NULL,
              multiple = TRUE
            ),
            shiny::checkboxInput(
              inputId = "toggle_all_schemes",
              label = bslib::tooltip(
                trigger = list(
                  "Select all schemes?",
                  bsicons::bs_icon("info-circle")
                ),
                "Select all available schemes at once."
              ),
              value = TRUE
            ),
            # add explanatory note
            shiny::markdown('##### Key'),
            shiny::markdown('✏️ = Scenario is not finalised')
          ),
          ### mitigator select ----
          bslib::accordion_panel(
            id = "accordion_mitigators",
            title = "Select mitigators",
            icon = bsicons::bs_icon("sliders"),
            shiny::selectInput(
              inputId = "activity_type",
              label = bslib::tooltip(
                trigger = list(
                  "Activity type",
                  bsicons::bs_icon("info-circle")
                ),
                "Select an activity type to filter the mitigators by."
              ),
              choices = NULL,
              selected = NULL,
              multiple = FALSE
            ),
            shiny::selectInput(
              inputId = "mitigator_type",
              label = bslib::tooltip(
                trigger = list(
                  "Mitigator type",
                  bsicons::bs_icon("info-circle")
                ),
                "Select the type of mitigator to filter by."
              ),
              choices = NULL,
              selected = NULL,
              multiple = FALSE
            ),
            shiny::selectInput(
              inputId = "mitigator_subset",
              label = bslib::tooltip(
                trigger = list(
                  "Mitigator subset",
                  bsicons::bs_icon("info-circle")
                ),
                "Select a subset group to filter mitigators by."
              ),
              choices = NULL,
              selected = NULL,
              multiple = FALSE
            ),
            shiny::selectInput(
              inputId = "mitigators",
              label = bslib::tooltip(
                trigger = list(
                  "Mitigators to visualise",
                  bsicons::bs_icon("info-circle")
                ),
                "Prepopulated given the mitigator filters applied above."
              ),
              choices = NULL,
              selected = NULL,
              multiple = TRUE
            )
          )
        )
      ),
      ## navpanel ----
      ### information ----
      bslib::nav_panel(
        id = "nav_panel_info",
        title = "Information",
        bslib::layout_column_wrap(
          bslib::card(
            id = "card_about",
            bslib::card_header("About"),
            md_file_to_html("app", "text", "about.md")
          ),
          bslib::card(
            id = "card_how_to_use",
            bslib::card_header("How to use"),
            md_file_to_html("app", "text", "how-to.md")
          )
        )
      ),
      ### point range plots ----
      bslib::nav_panel(
        id = "nav_panel_pointrange",
        title = "Point-ranges",
        bslib::card(
          full_screen = TRUE,
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              title = "Point-range settings",
              open = TRUE,
              width = 260,
              shiny::checkboxInput(
                inputId = "toggle_horizon_pointrange",
                label = bslib::tooltip(
                  trigger = list(
                    "Standardise by horizon length?",
                    bsicons::bs_icon("info-circle")
                  ),
                  "Divides the scheme's chosen mitigator values by the number of years between the chosen start and final year."
                ),
                value = FALSE
              ),
              shiny::checkboxInput(
                inputId = "toggle_mitigator_code_pointrange",
                label = bslib::tooltip(
                  trigger = list(
                    "Show mitigator code?",
                    bsicons::bs_icon("info-circle")
                  ),
                  "Replaces the full mitigator name with the mitigator code."
                ),
                value = FALSE
              ),
              shiny::checkboxInput(
                inputId = "toggle_invert_facets",
                label = bslib::tooltip(
                  trigger = list(
                    "Facet by scheme?",
                    bsicons::bs_icon("info-circle")
                  ),
                  "Invert the pointrange plots to show mitigators on the y axis and scheme as the faceting variable."
                ),
                value = FALSE
              ),
              shiny::checkboxInput(
                inputId = "toggle_nee_reference_range",
                label = bslib::tooltip(
                  trigger = list(
                    "Show NEE range?",
                    bsicons::bs_icon("info-circle")
                  ),
                  "Include reference results from the National Elicitation Exercise (NEE). These values are shown as horizontal bars behind each point illustrating the 10% to 90% range, with a vertical line marking the mean value."
                ),
                value = TRUE
              ),
              shiny::sliderInput(
                inputId = "facet_columns",
                label = bslib::tooltip(
                  trigger = list(
                    "Number of facet columns",
                    bsicons::bs_icon("info-circle")
                  ),
                  "Choose the number of columns over which to break the faceted pointrange charts."
                ),
                min = 1,
                max = 5, # will be reactively updated to match the number of facets
                step = 1,
                value = 5,
                round = TRUE,
                ticks = FALSE
              ),
              shiny::bookmarkButton(
                label = "Bookmark ",
                icon = shiny::icon("bookmark", lib = "glyphicon"),
                style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
              )
            ),
            shiny::plotOutput("pointrange"),
          )
        )
      ),
      ### heatmap plots ----
      bslib::nav_panel(
        id = "nav_panel_heatmaps",
        title = "Heatmaps",

        bslib::navset_card_underline(
          id = 'nav_panel_heatmaps_tabs',
          full_screen = TRUE,

          #### heatmaps ----
          bslib::nav_panel(
            title = bslib::tooltip(
              trigger = list(
                'Heatmaps',
                bsicons::bs_icon('info-circle')
              ),
              'Customisable heatmaps showing distributions of values by mitigator and scheme'
            ),
            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                title = "Heatmap settings",
                open = TRUE,
                shiny::selectInput(
                  inputId = "heatmap_type",
                  label = bslib::tooltip(
                    trigger = list(
                      "Value type",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Schemes' low or high 80% confidence internal selection in the NHP inputs app, or the range or midpoint of these."
                  ),
                  choices = c(
                    Binary = "value_binary",
                    Midpoint = "value_mid",
                    Range = "value_range",
                    Low = "value_lo",
                    High = "value_hi"
                  ),
                  selected = "value_mid",
                  multiple = FALSE
                ),
                shiny::checkboxInput(
                  inputId = "toggle_horizon_heatmap",
                  label = bslib::tooltip(
                    trigger = list(
                      "Standardise by horizon length?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Divides the scheme's chosen mitigator values by the number of years between the chosen start and final year."
                  ),
                  value = FALSE
                ),
                shiny::checkboxInput(
                  inputId = "toggle_mitigator_name",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show mitigator names?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Plots mitigator names on the y-axis (default) or switch off to display mitigator codes instead."
                  ),
                  value = TRUE
                ),
                shiny::bookmarkButton(
                  label = "Bookmark",
                  icon = shiny::icon("bookmark", lib = "glyphicon"),
                  style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                )
              ),
              shiny::plotOutput("heatmap"),
            )
          ),

          #### mitigator coverage ----
          bslib::nav_panel(
            title = bslib::tooltip(
              trigger = list(
                'Mitigator coverage',
                bsicons::bs_icon('info-circle')
              ),
              'The proportion of schemes using each mitigator',
            ),
            DT::DTOutput("mitigator_uptake_dt")
          ),

          #### scheme coverage ----
          bslib::nav_panel(
            title = bslib::tooltip(
              trigger = list(
                'Scheme coverage',
                bsicons::bs_icon('info-circle')
              ),
              'The proportion of mitigators in use by each scheme. Selected schemes are shown in bold, the focal scheme is highlighted in red.',
            ),
            DT::DTOutput("scheme_uptake_dt")
          ),
        ),
      ),

      ### data -----
      bslib::nav_panel(
        id = "nav_panel_data",
        title = "Data",
        bslib::navset_card_underline(
          full_screen = TRUE,
          bslib::nav_panel(
            "Raw data",
            DT::DTOutput("raw_data_dt")
          ),
          bslib::nav_panel(
            "Mitigator lookup",
            DT::DTOutput("mitigator_lookup_dt")
          ),
          bslib::nav_panel(
            "Scheme lookup",
            DT::DTOutput("scheme_lookup_dt")
          )
        )
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "NHP mitigator comparisons"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
