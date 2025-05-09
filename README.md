# nhp_inputs_report_app

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

## Purpose

An app built with [{shiny}](https://shiny.posit.co/), [{golem}](https://thinkr-open.github.io/golem/) and [{bslib}](https://rstudio.github.io/bslib/) to explore and compare schemes' mitigator selections as part of the National Hospital Programme (NHP) modelling process. 

This tool is designed primarily for use by model-relationship managers (MRMs) in discussion with schemes so that the mitigator selections can be refined before being finalised.

The app is [deployed to Posit Connect](https://connect.strategyunitwm.nhs.uk/nhp/mitigator-comparisons/) (login/permissions required).

## Run the app

### Deploy

You can redeploy the app to Posit Connect using the `dev/03_deploy.R` script.
This usually happens after a new GitHub release/Git tag.

The script checks for the 'app ID' in the `rsconnect/` folder of your local project root, which is generated when you first deploy.
Otherwise you can find the ID by opening the app from the Posit Connect 'Content' page and then looking for 'Content ID' in the Settings > Info panel of the interface.

### Locally

You can run the app locally, but some environmental variables are needed to fetch data, for example.
You will need to add a `.Renviron` file to the project directory that contains the variables named in the `.Renviron.example` file.
You can ask a member of the Data Science team for the values to populate this file.
Remember to restart your session after you've updated your `.Renviron` file.

## Data

### Mitigators

This app fetches model results files from Azure. 
These files are large json files that bundle the model results and the the mitigator selections, which are stored in an element called 'params' (parameters).

To avoid the app having to read the entire json file for each scheme's model scenario, there is a system to pre-prepare the params alone.
[A scheduled Quarto document on Posit Connect](https://connect.strategyunitwm.nhs.uk/nhp/tagged-runs-params-report/) has code to select the appropriate json file, extract the params and [save them as an RDS file to a pin](https://connect.strategyunitwm.nhs.uk/content/32c7f642-e420-448d-b888-bf655fc8fa8b/) on Posit Connect.
It also saves [a CSV file to another pin](https://connect.strategyunitwm.nhs.uk/content/811dbaf9-18fe-43aa-bf8e-06b0df66004e/) that contains metadata about the model scenarios.
The app then reads data from these pins using [the {pins} package](https://pins.rstudio.com/).
A [blogpost by the Data Science team](https://the-strategy-unit.github.io/data_science/blogs/posts/2024-05-22-storing-data-safely/#posit-connect-pins) contains note on authenticating RStudio with Posit Connect, should you need to.

Schemes run many scenarios, but the app only displays data from a single scenario, preferably the one used to compile their outputs ('final') report.
MRMs tell the Data Science team which particular results file should be labelled manually on Azure with the 'run stage' metadata of 'final' (and possibly 'intermediate' or 'initial').
There's [a handy lookup table](https://connect.strategyunitwm.nhs.uk/nhp/tagged_runs/nhp-tagged-runs.html) (login/permissions required) where you can see the scenario files that have been given run-stage metadata.

### Supporting

Supporting data is fetched from Azure.
This includes lookups for mitigators and schemes, baseline trend data, as well as data from the National Elicitation Exercise (NEE).
