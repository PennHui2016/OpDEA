library(shinydashboard)
library(shiny)
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
source('./R/run_DEA_single.R')
source('./R/preprocessing_pro_counts.R')
source('./R/preprocessing_pro_intensity.R')
source('./R/preprocessing_pro_intensity_DEqMS.R')
source('./R/preprocessing_pro_intensity_DEqMS.R')
source('./R/run_DEA_single.R')
source('./R/run_DEA_ens.R')
source('./R/iq-fast_MOD.R')
library(iq)
source('./R/limma_DE.R')
source('./R/ROTS_DE.R')
source('./R/DEP_DE.R')
source('./R/proDA_DE.R')
source('./R/DEqMS_DE.R')
source('./R/plgem_DE.R')
source('./R/edgeR_DE.R')
source('./R/beta_binomial_DE.R')
source('./R/ANOVA_DE.R')
source('./R/SAM_DE.R')
source('./R/ttest_DE.R')
source('./R/MSstats_DE.R')
source('./R/ensemble.R')

#source('./R/app_config.R')
#addResourcePath('img', 'app/img')
addResourcePath('img', 'inst/app/img')
#addResourcePath('www', 'app/www')
#addResourcePath('root_fold', 'inst/app/root_fold')

add_resource_path(
  "www",
  app_sys("app/www")
)

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
