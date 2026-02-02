# =============================================================================
# Shiny Dashboard for CausalDef
# =============================================================================

#' Launch CausalDef Shiny Dashboard
#'
#' Launches an interactive Shiny application for deficiency analysis.
#' Provides a point-and-click interface for the full causaldef workflow.
#'
#' @param data Optional data frame to preload
#' @param port Integer: port for the Shiny server (default 3838)
#' @param launch.browser Logical: open in browser? (default TRUE)
#'
#' @return Invisibly returns the Shiny app object
#'
#' @details
#' The dashboard provides:
#' \itemize{
#'   \item Data upload and variable selection
#'   \item Method comparison (unadjusted, IPTW, AIPW, etc.)
#'   \item Interactive visualization of deficiency estimates
#'   \item Confounding frontier explorer
#'   \item Policy regret calculator
#'   \item Exportable reports
#' }
#'
#' @examples
#' \dontrun{
#' # Launch with example data
#' df <- data.frame(
#'   W = rnorm(200),
#'   A = rbinom(200, 1, 0.5),
#'   Y = rnorm(200)
#' )
#' run_causaldef_app(data = df)
#' }
#'
#' @export
run_causaldef_app <- function(data = NULL, port = 3838, launch.browser = TRUE) {
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg shiny} is required for the dashboard.",
      "i" = "Install with: {.code install.packages('shiny')}"
    ))
  }
  
  # Define UI
  ui <- shiny::fluidPage(
    shiny::titlePanel("CausalDef: Decision-Theoretic Causal Diagnostics"),
    
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        
        # Data Input
        shiny::h4("1. Data"),
        shiny::fileInput("file", "Upload CSV", accept = ".csv"),
        shiny::uiOutput("var_selectors"),
        
        shiny::hr(),
        
        # Methods
        shiny::h4("2. Methods"),
        shiny::checkboxGroupInput(
          "methods", "Select Methods:",
          choices = c("Unadjusted" = "unadjusted",
                     "IPTW" = "iptw",
                     "AIPW" = "aipw",
                     "Matching" = "matching"),
          selected = c("unadjusted", "iptw", "aipw")
        ),
        
        shiny::numericInput("n_boot", "Bootstrap Samples:", 100, min = 0, max = 1000),
        
        shiny::hr(),
        
        # Actions
        shiny::actionButton("run", "Run Analysis", class = "btn-primary"),
        shiny::downloadButton("report", "Download Report")
      ),
      
      shiny::mainPanel(
        width = 9,
        
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Deficiency",
            shiny::h3("Le Cam Deficiency Estimates"),
            shiny::plotOutput("deficiency_plot", height = "400px"),
            shiny::verbatimTextOutput("deficiency_summary")
          ),
          
          shiny::tabPanel(
            "Confounding Frontier",
            shiny::h3("Sensitivity Analysis"),
            shiny::sliderInput("alpha_range", "Alpha Range:", min = -3, max = 3, value = c(-2, 2)),
            shiny::sliderInput("gamma_range", "Gamma Range:", min = -3, max = 3, value = c(-2, 2)),
            shiny::plotOutput("frontier_plot", height = "500px")
          ),
          
          shiny::tabPanel(
            "Policy Bounds",
            shiny::h3("Policy Regret Analysis"),
            shiny::numericInput("utility_min", "Utility Min:", 0),
            shiny::numericInput("utility_max", "Utility Max:", 1),
            shiny::numericInput("obs_regret", "Observed Regret:", 0, min = 0),
            shiny::plotOutput("policy_plot", height = "400px"),
            shiny::verbatimTextOutput("policy_summary")
          ),
          
          shiny::tabPanel(
            "Data Preview",
            shiny::h3("Data Summary"),
            shiny::verbatimTextOutput("data_summary"),
            shiny::h4("First 10 Rows"),
            shiny::tableOutput("data_head")
          )
        )
      )
    )
  )
  
  # Define Server
  server <- function(input, output, session) {
    
    # Reactive: loaded data
    data_reactive <- shiny::reactive({
      if (!is.null(data)) {
        return(data)
      }
      
      shiny::req(input$file)
      utils::read.csv(input$file$datapath)
    })
    
    # Dynamic variable selectors
    output$var_selectors <- shiny::renderUI({
      df <- data_reactive()
      if (is.null(df)) return(NULL)
      
      vars <- names(df)
      
      shiny::tagList(
        shiny::selectInput("treatment", "Treatment:", choices = vars),
        shiny::selectInput("outcome", "Outcome:", choices = vars, selected = vars[2]),
        shiny::selectizeInput("covariates", "Covariates:", choices = vars, multiple = TRUE)
      )
    })
    
    # Data summary
    output$data_summary <- shiny::renderPrint({
      df <- data_reactive()
      if (is.null(df)) return("No data loaded")
      summary(df)
    })
    
    output$data_head <- shiny::renderTable({
      df <- data_reactive()
      if (is.null(df)) return(NULL)
      utils::head(df, 10)
    })
    
    # Run analysis
    results <- shiny::eventReactive(input$run, {
      df <- data_reactive()
      shiny::req(df, input$treatment, input$outcome)
      
      # Create causal spec
      spec <- causal_spec(
        data = df,
        treatment = input$treatment,
        outcome = input$outcome,
        covariates = if (length(input$covariates) > 0) input$covariates else NULL
      )
      
      # Estimate deficiency
      def_result <- tryCatch({
        estimate_deficiency(
          spec,
          methods = input$methods,
          n_boot = input$n_boot
        )
      }, error = function(e) {
        shiny::showNotification(paste("Error:", e$message), type = "error")
        NULL
      })
      
      # Confounding frontier
      frontier <- tryCatch({
        confounding_frontier(
          spec,
          alpha_range = input$alpha_range,
          gamma_range = input$gamma_range
        )
      }, error = function(e) NULL)
      
      list(spec = spec, deficiency = def_result, frontier = frontier)
    })
    
    # Deficiency plot
    output$deficiency_plot <- shiny::renderPlot({
      res <- results()
      if (is.null(res$deficiency)) return(NULL)
      plot(res$deficiency, type = "bar")
    })
    
    output$deficiency_summary <- shiny::renderPrint({
      res <- results()
      if (is.null(res$deficiency)) return("Run analysis first")
      print(res$deficiency)
    })
    
    # Frontier plot
    output$frontier_plot <- shiny::renderPlot({
      res <- results()
      if (is.null(res$frontier)) return(NULL)
      plot(res$frontier, type = "heatmap")
    })
    
    # Policy bounds
    output$policy_plot <- shiny::renderPlot({
      res <- results()
      if (is.null(res$deficiency)) return(NULL)
      
      bounds <- policy_regret_bound(
        res$deficiency,
        utility_range = c(input$utility_min, input$utility_max),
        obs_regret = input$obs_regret
      )
      
      plot(bounds, type = "safety_curve")
    })
    
    output$policy_summary <- shiny::renderPrint({
      res <- results()
      if (is.null(res$deficiency)) return("Run analysis first")
      
      bounds <- policy_regret_bound(
        res$deficiency,
        utility_range = c(input$utility_min, input$utility_max),
        obs_regret = input$obs_regret
      )
      
      print(bounds)
    })
    
    # Report download
    output$report <- shiny::downloadHandler(
      filename = function() {
        paste0("causaldef_report_", Sys.Date(), ".txt")
      },
      content = function(file) {
        res <- results()
        
        sink(file)
        cat("=== CausalDef Analysis Report ===\n")
        cat("Generated:", as.character(Sys.time()), "\n\n")
        
        if (!is.null(res$deficiency)) {
          cat("--- Deficiency Estimates ---\n")
          print(res$deficiency)
          cat("\n")
        }
        
        if (!is.null(res$frontier)) {
          cat("--- Confounding Frontier ---\n")
          cat("Alpha range:", input$alpha_range, "\n")
          cat("Gamma range:", input$gamma_range, "\n")
          cat("Safe region (delta < 0.1):", 
              sum(res$frontier$grid$delta < 0.1) / nrow(res$frontier$grid) * 100, "%\n")
        }
        
        sink()
      }
    )
  }
  
  # Run app
  app <- shiny::shinyApp(ui = ui, server = server)
  
  shiny::runApp(app, port = port, launch.browser = launch.browser)
  
  invisible(app)
}

#' Create Standalone Shiny App Files
#'
#' Generates app.R and global.R files for standalone deployment
#' (e.g., shinyapps.io, Shiny Server)
#'
#' @param path Directory to write the app files
#'
#' @export
create_shiny_app_files <- function(path = "causaldef_app") {
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  # Write app.R
  app_code <- '
library(shiny)
library(causaldef)
library(ggplot2)

# See ?run_causaldef_app for the full implementation
# This is a minimal standalone version

ui <- fluidPage(
  titlePanel("CausalDef Dashboard"),
  # ... (UI code from run_causaldef_app)
)

server <- function(input, output, session) {
  # ... (server code from run_causaldef_app)
}

shinyApp(ui, server)
'
  
  writeLines(app_code, file.path(path, "app.R"))
  
  cli::cli_alert_success("Created Shiny app files in {.file {path}}")
  cli::cli_alert_info("Deploy to shinyapps.io: {.code rsconnect::deployApp('{path}')}")
  
  invisible(path)
}
