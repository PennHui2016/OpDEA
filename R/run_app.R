#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options

library(golem)
source("./R/plot_figures.R")
source('./R/tabs_contents.R')
source("./R/Data_access.R")
source("./R/helps.R")
source("./R/suggestions.R")
source('./R/response_suggest.R')
source('./R/app_server.R')
source('./R/app_ui.R')
#source('./R/app_config.R')

addResourcePath('img', 'inst/app/img')
#addResourcePath('root_fold', 'inst/app/root_fold')
run_app <- function(
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  ...
) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
