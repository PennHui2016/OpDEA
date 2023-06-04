#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
library(shinydashboard)
library(shiny)
library(threejs)
library(DT)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggsci)
library("readxl")
library(ggalluvial)
app_ui <- function(request) {
  dashboardPage(
    dashboardHeader(title = "OpDEA"),
    dashboardSidebar(
      #uiOutput("sidebarControls"),
      sidebarMenuOutput("menu")
    ),
    dashboardBody(
      tabItems(
        tabItem("front",
                h3("OpDEA: Achieving Superior Differential Analysis in Proteomics Workflows via Ensemble Inference Optimization"),
                br(),
                tabIntro,
                tabWf_intro,
                tabLOPOCV,
                tabcls,
                tabAC
        ),
        tabItem("models1",
                h3("FragPipe workflow benchmarking results"),
                br(),
                tab1
        ),
        tabItem("models2",
                h3("maxquant workflow benchmarking results"),
                br(),
                tab2
        ),
        tabItem("models3",
                h3("DIA-NN workflow benchmarking results"),
                br(),
                tab3
        ),
        tabItem("sug1",
                h3("Recommend optimal Workflow for DDA label-free data"),
                br(),
                tabpanel1,

        ),
        tabItem("sug2",
                h3("Recommend optimal Workflow for DIA label-free data"),
                br(),
                tabpanel2
        ),
        tabItem("accData",
                h3("Avaliable Data"),
                br(),
                Datas1,#raw data link
                Datas2,#quantification results
                Datas3#our benchmarking results
                #Data4#ranks of workflows

        ),
        tabItem("Help1",
                h3("View benchmarking results"),
                br(),
                help_tab1
        ),
        tabItem("Help2",
                h3("Optimal workflow recommendation"),
                br(),
                help_tab2
        ),
        tabItem("Help3",
                h3("Contact Info"),
                br(),
                help_tab3

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
  #addResourcePath("inst/app/root_fold/workflows.xlsx", "inst/app/root_fold/workflows.xlsx")
  addResourcePath('img', 'inst/app/img')
  addResourcePath('root_fold', 'inst/app/root_fold')
  add_resource_path(
    "root_fold",
    app_sys("app/root_fold")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "ShinyDektop"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
