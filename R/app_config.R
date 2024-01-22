#' Access files in the current app
#'
#' NOTE: If you manually change your package name in the DESCRIPTION,
#' don't forget to change it here too, and in the config file.
#' For a safer name change mechanism, use the `golem::set_golem_name()` function.
#'
#' @param ... character vectors, specifying subdirectory and file(s)
#' within your package. The default, none, returns the root of the app.
#'
#' @noRd
app_sys <- function(...) {
  system.file(..., package = "OpDEA")
}


#diann_rk<-app_sys('diann_spr_avg_rank_test.csv')
#fg_rk<-app_sys('fragpipe_spr_avg_rank_test.csv')
#mq_rk<-app_sys('maxquant_spr_avg_rank_test.csv')
#me_diann<-app_sys('metrics_diann.xlsx')
#me_fg<-app_sys('metrics_fragpipe.xlsx')
#me_max<-app_sys('metrics_maxquant.xlsx')
#sum_cat<-app_sys('sum_CatBoost.xlsx')
#' Read App Config
#'
#' @param value Value to retrieve from the config file.
#' @param config GOLEM_CONFIG_ACTIVE value. If unset, R_CONFIG_ACTIVE.
#' If unset, "default".
#' @param use_parent Logical, scan the parent directory for config file.
#' @param file Location of the config file
#'
#' @noRd
get_golem_config <- function(
  value,
  config = Sys.getenv(
    "GOLEM_CONFIG_ACTIVE",
    Sys.getenv(
      "R_CONFIG_ACTIVE",
      "default"
    )
  ),
  use_parent = TRUE,
  # Modify this if your config file is somewhere else
  file = app_sys("golem-config.yml")
) {
  config::get(
    value = value,
    config = config,
    file = file,
    use_parent = use_parent
  )
}
