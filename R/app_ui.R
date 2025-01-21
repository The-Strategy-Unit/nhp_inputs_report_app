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
    shinybusy::add_busy_spinner(
      spin = "double-bounce",
      margins = c(10, 10)
    ),
    # Application UI logic
    bslib::page_navbar(
      theme = bslib::bs_theme(version = 5L),
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
                "The scheme to highlight in plots. Causes autoselection of peers."
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
                "Automatically select all schemes at once."
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

            # new mitigator selection ---
            shinyWidgets::panel(
              # NB, using functions from the github {datamods} package repo
              # as using the CRAN version results in certain bslib containers to
              # lose page-filling ability. The new functions in the github
              # repo do not do this. It is possible the CRAN version will catch-up
              # at some point so we can remove the duplicated functions.
              # NB, datamod functions are duplicated in the following files:
              # datamods_select_group.R
              # datamods_utils_shiny.R
              # datamods_utils.R
              select_group_ui(
                id = 'mitigator_filters',
                vs_args = list(
                  search = TRUE,
                  showValuesAsTags = TRUE
                ),
                inline = FALSE,
                params = list(
                  mitigator_type = list(
                    inputId = 'mitigator_type',
                    label = 'Mitigator type',
                    placeholder = '...'
                  ),
                  activity_type = list(
                    inputId = 'activity_type',
                    label = 'Activity type',
                    placeholder = '...'
                  ),
                  grouping = list(
                    inputId = 'grouping',
                    label = 'Mitigator group',
                    placeholder = '...'
                  ),
                  strategy_subset = list(
                    inputId = 'strategy_subset',
                    label = 'Strategy subset',
                    placeholder = '...'
                  ),
                  mitigator_name = list(
                    inputId = 'mitigator_name',
                    label = 'Mitigator',
                    placeholder = '...'
                  )
                )
              ),

              # add some separation between 'reset filters' and below button
              shiny::br(),
              shiny::br(),

              # button to add mitigators to the selection box (🍫)
              shiny::actionButton(
                inputId = 'mitigators_add_to_selected',
                label = 'Add to selection',
                icon = shiny::icon('arrow-down', lib = 'font-awesome'),
                class = 'btn-light'
              )
            ),

            # list mitigators
            shiny::selectizeInput(
              inputId = 'mitigators',
              label = 'Selected mitigators',
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list('plugins' = list('remove_button'))
            ),

            # clear selected mitigators
            shiny::actionLink(
              inputId = 'clear_selected_mitigators',
              label = bslib::tooltip(
                trigger = list(
                  bsicons::bs_icon("x-lg"),
                  "Clear selected"
                ),
                "Remove all mitigators from the current selection."
              ),
            ),
          ),
          ### other settings ----
          bslib::accordion_panel(
            id = "accordion_other_settings",
            title = "Other settings",
            icon = bsicons::bs_icon("gear"),
            shinyWidgets::radioGroupButtons(
              inputId = "values_displayed",
              label = bslib::tooltip(
                trigger = list(
                  "Values to display",
                  bsicons::bs_icon("info-circle")
                ),
                "Select whether values are shown as the amount mitigated or the expected activity levels following mitigation."
              ),
              choices = c("Percent of activity mitigated", "80% prediction interval"),
              selected = "Percent of activity mitigated",
              direction = "vertical",
              checkIcon = list(
                yes = icon("ok", lib = "glyphicon")
              )
            ),
            shinyWidgets::materialSwitch(
              inputId = "standardise_2041",
              label = bslib::tooltip(
                trigger = list(
                  "Standardise to 2041?",
                  bsicons::bs_icon("info-circle")
                ),
                "Standardise values by extrapolating linearly to 2041, which makes a direct comparison easier."
              ),
              value = FALSE,
              status = "primary",
              right = TRUE
            ),
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

        bslib::navset_card_underline(
          id = 'nav_panel_heatmaps_tabs',
          full_screen = TRUE,

          #### pointrange ----
          bslib::nav_panel(
            title = bslib::tooltip(
              trigger = list(
                'Point-range',
                bsicons::bs_icon('info-circle')
              ),
              'Customisable point-ranges showing distributions of values by mitigator and scheme'
            ),
            #fillable = TRUE,
            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                #title = "Point-range settings",
                open = TRUE,
                width = 350,
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
                    "Include reference results from the National Elicitation Exercise (NEE) for 2039/40. These values are shown as horizontal bars behind each point illustrating the 10% to 90% range, with a vertical line marking the mean value."
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
                shinyWidgets::materialSwitch(
                  inputId = "toggle_aggregate_summary",
                  label = bslib::tooltip(
                    trigger = list(
                      "Summary",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Add an aggregate summary of selected schemes' responses to the plot (includes all but the focal scheme)."
                  ),
                  value = FALSE,
                  status = "primary",
                  right = TRUE
                ),
                shiny::checkboxInput(
                  inputId = "toggle_aggregate_summary_minmaxrange",
                  label = bslib::tooltip(
                    trigger = list(
                      "Summary full range",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Show the full range of selected schemes' responses, i.e. the extreme upper and lower values, or switch off to view the average (mean) range."
                  ),
                  value = FALSE
                ),
                shiny::bookmarkButton(
                  label = "Bookmark ",
                  icon = shiny::icon("bookmark", lib = "glyphicon"),
                  style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                ),

                bslib::accordion(
                  open = FALSE,
                  bslib::accordion_panel(
                    title = "About",
                    md_file_to_html("app", "text", "about_pointrange.md")
                  )
                )
              ), # end sidebar
              shiny::plotOutput("pointrange")
            )
          ),

          #### densities -----
          bslib::nav_panel(
            title = bslib::tooltip(
              trigger = list(
                'Distributions',
                bsicons::bs_icon('info-circle')
              ),
              'Mixture distributions of values provided by schemes for mitigators'
            ),

            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                #title = "Density settings",
                open = TRUE,
                width = 350,

                shiny::checkboxInput(
                  inputId = "toggle_mixture_distribution_ecdf",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show ECDF?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Plot the Empirical Cumulative Distribution Function (ECDF), or leave unchecked for the Probability Density Function (PDF)."
                  ),
                  value = FALSE
                ),

                shiny::checkboxInput(
                  inputId = "toggle_nee_reference_range_density",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show NEE range?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Include reference results from the National Elicitation Exercise (NEE) for 2039/40. These values are shown as horizontal bars behind each point illustrating the 10% to 90% range, with a vertical line marking the mean value."
                  ),
                  value = TRUE
                ),

                shiny::bookmarkButton(
                  label = "Bookmark ",
                  icon = shiny::icon("bookmark", lib = "glyphicon"),
                  style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                ),

                bslib::accordion(
                  open = FALSE,
                  bslib::accordion_panel(
                    title = "About",
                    md_file_to_html("app", "text", "about_distributions.md")
                  )
                )
              ), # end sidebar
              shiny::plotOutput('mixture_distributions')
            )
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
                width = 350,
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
                ),
                bslib::accordion(
                  open = FALSE,
                  bslib::accordion_panel(
                    title = "About",
                    md_file_to_html("app", "text", "about_heatmaps.md")
                  )
                )
              ), # end sidebar
              shiny::plotOutput("heatmap")
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
            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                open = TRUE,
                width = 350,

                bslib::accordion(
                  open = FALSE,
                  bslib::accordion_panel(
                    title = "About",
                    md_file_to_html("app", "text", "about_mitigator_coverage.md")
                  )
                )
              ), # end sidebar
              DT::DTOutput("mitigator_uptake_dt"),
            ),
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
            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                width = 350,
                open = TRUE,

                bslib::accordion(
                  open = FALSE,
                  bslib::accordion_panel(
                    title = "About",
                    md_file_to_html("app", "text", "about_scheme_coverage.md")
                  )
                )
              ), # end sidebar
              DT::DTOutput("scheme_uptake_dt")
            ),
          ),
        ),
      ),

      ### contextual plots ----
      bslib::nav_panel(
        id = "nav_panel_contextual",
        title = "Contextual",

        bslib::navset_card_underline(
          full_screen = TRUE,

          #### baseline plot ----
          bslib::nav_panel(
            title = bslib::tooltip(
              trigger = list(
                "Baseline comparison",
                bsicons::bs_icon("info-circle")
              ),
              "A scatter plot comparing schemes' baseline values on the x-axis and their mitigator value inputs on the y-axis."
            ),
            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                width = 350,
                open = TRUE,

                shiny::checkboxInput(
                  inputId = "toggle_contextual_baseline_range",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show 80% range?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Shows the 80% range between upper and lower values as points connected by a line, disable to view just the mid-points."
                  ),
                  value = TRUE
                ),
                shiny::checkboxInput(
                  inputId = "toggle_contextual_baseline_quadrants",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show quadrant lines?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Shows dotted lines indicating schemes' average (mean) values for each axis, dividing the plot into quadrants."
                  ),
                  value = TRUE
                ),
                shiny::checkboxInput(
                  inputId = "toggle_contextual_baseline_schemecode",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show scheme codes?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Shows a label for each scheme code on the plot."
                  ),
                  value = FALSE
                ),
                shiny::sliderInput(
                  inputId = "slider_contextual_baseline_height",
                  label = bslib::tooltip(
                    trigger = list(
                      "Plot height",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Set the height of each plot in pixels."
                  ),
                  value = 250,
                  min = 150,
                  max = 400,
                  step = 50
                ),
                shiny::bookmarkButton(
                  label = "Bookmark",
                  icon = shiny::icon("bookmark", lib = "glyphicon"),
                  style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                ),
                bslib::accordion(
                  open = FALSE,
                  bslib::accordion_panel(
                    title = "About",
                    md_file_to_html("app", "text", "about_baseline.md")
                  )
                )
              ), # end sidebar
              plotly::plotlyOutput("contextual_baseline")
            )
          ),

          #### trendline plot ----
          bslib::nav_panel(
            title = bslib::tooltip(
              trigger = list(
                "Trendline comparison",
                bsicons::bs_icon("info-circle")
              ),
              "A timeseries plot showing historical rates of activity for the focal scheme enabling comparisons with their predicted activity."
            ),
            bslib::layout_sidebar(
              sidebar = bslib::sidebar(
                width = 350,
                open = TRUE,

                shiny::checkboxInput(
                  inputId = "toggle_contextual_trendline_otherschemes",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show other schemes?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Shows trendline plots for schemes other than the focal scheme as grey traces."
                  ),
                  value = FALSE
                ),
                shiny::checkboxInput(
                  inputId = "toggle_contextual_trendline_horizon_timeline",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show horizon on timeline?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Shows the predicted activity at the horizon year on the timeline. The predicted interval is shown as a point range on the horizon year with dotted lines connecting it with the baseline rate."
                  ),
                  value = TRUE
                ),
                shiny::checkboxInput(
                  inputId = "toggle_contextual_trendline_horizon_overlay",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show horizon as overlay?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Shows the predicted activity as an overlay. The predicted interval is shown as three lines (low, mid and high) with a coloured ribbon across the historical activity timeline."
                  ),
                  value = FALSE
                ),
                shiny::checkboxInput(
                  inputId = "toggle_contextual_trendline_average",
                  label = bslib::tooltip(
                    trigger = list(
                      "Show pre-baseline average?",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Shows the average (mean) rate up to the baseline year along with a window of two standard deviations above and below to indicate the pre-baseline activity and range."
                  ),
                  value = FALSE
                ),
                shiny::sliderInput(
                  inputId = "slider_contextual_trendline_height",
                  label = bslib::tooltip(
                    trigger = list(
                      "Plot height",
                      bsicons::bs_icon("info-circle")
                    ),
                    "Set the height of each plot in pixels."
                  ),
                  value = 400,
                  min = 250,
                  max = 600,
                  step = 50
                ),
                shiny::bookmarkButton(
                  label = "Bookmark",
                  icon = shiny::icon("bookmark", lib = "glyphicon"),
                  style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                ),
                bslib::accordion(
                  open = FALSE,
                  bslib::accordion_panel(
                    title = "About",
                    md_file_to_html("app", "text", "about_trendline.md")
                  )
                )
              ), # end sidebar
              plotly::plotlyOutput("contextual_trendline")
            )
          )
        )
      ),

      ### data -----
      bslib::nav_panel(
        id = "nav_panel_data",
        title = "Data",
        bslib::navset_card_underline(
          full_screen = TRUE,
          bslib::nav_panel(
            title = bslib::tooltip(
              trigger = list(
                "Raw data",
                bsicons::bs_icon('info-circle')
              ),
              "An interactive table of underlying data. Contains scheme, mitigator and model-run metadata; schemes' selected mitigator values; and results of the National Elicitation Exercise (NEE). Column-name prefixes are 'pm' for percent mitigated and 'pi' for 80% prediction interval. Suffixes include 'p10' to mean the 10th percentile.",
            ),
            DT::DTOutput("raw_data_dt")
          ),
          bslib::nav_panel(
            title = "Mitigator lookup",
            DT::DTOutput("mitigator_lookup_dt")
          ),
          bslib::nav_panel(
            title = "Scheme lookup",
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
