# ====================================================================
# prepare
# ====================================================================
library(shiny)
library(shinydashboard)
library(DT)
library(dplyr)
library(shinyWidgets)

# Get list of available proteins
# protein_files <- list.files("data/1aa_del", pattern = "\\.csv$", full.names = TRUE)
protein_files <- list.files("/srv/shiny-server/db/pon_del/data/1aa_del/output", pattern = "\\.csv$", full.names = TRUE)
protein_names <- basename(protein_files)
protein_display_names <- sub("\\.csv$", "", protein_names)

# Load ID mapping dictionary
id_mapping <- read.csv("/srv/shiny-server/db/mapping/mane_dict.csv", stringsAsFactors = FALSE)
id_mapping = id_mapping%>%select(-NC_id)%>%filter(NP_id %in% protein_display_names)

# Demo inputs
# protein_demo <- "NP_id,start,end\nNP_000008.1,104,104\nNP_000009.1,299,299\nNP_000009.1,276,283\nNP_000009.1,504,505\nNP_570602.2,3,3\nNP_570602.2,4,5"
# transcript_demo <- "ENST_id,start,end\nENST00000242592.9,310,312\nENST00000356839.10,896,898\nENST00000356839.10,826,849\nENST00000356839.10,1511,1516\nENST00000263100.8,7,9\nENST00000263100.8,10,15"
# genomic_demo <- "chr,start,end\n12,120737085,120737087\n17,7222684,7222686\n17,7222250,7222273\n17,7224222,7224227\n19,58353429,58353431\n19,58353425,58353430"

protein_demo <- "NP_id,start,end\nNP_000008.1,104,104\nNP_000009.1,299,299"
transcript_demo <- "ENST_id,start,end\nENST00000356839.10,826,849\nENST00000356839.10,1511,1516"
genomic_demo <- "chr,start,end\n19,58353429,58353431\n19,58353425,58353430"

# Create a list of available ID types
id_types <- c(
  "Gene symbol" = "gene_symbol",
  "RefSeq protein" = "NP_id",
  "RefSeq transcript" = "NM_id",
  "Ensembl protein" = "ENSP_id",
  "Ensembl transcript" = "ENST_id",
  "Ensembl gene" = "ENSG_id"
)
  

ui = shinyUI(navbarPage(
  title = "PON-Del",
  tags$head(
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css"),
    tags$link(rel = "icon", type = "image/x-icon", href = "data:image/svg+xml,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><text y='.9em' font-size='90'>🧬</text></svg>"),
    tags$style(HTML("
      /* Navbar icon removed */
    "))
  ),
  tags$style(HTML("
    #protein_plot img {
      width: 100% !important;
      max-width: 100% !important;
      height: auto !important;
      display: block;
      margin-left: auto;
      margin-right: auto;
    }
    .well {
      width: 100% !important;
      max-width: 100% !important;
    }
    #protein_table table.dataTable {
      font-size: 13px !important;
    }
    #prediction_log {
      font-family: monospace;
      background-color: #f8f9fa;
      padding: 10px;
      border-radius: 4px;
      max-height: 200px;
      overflow-y: auto;
      white-space: pre-wrap;
      font-size: 12px;
    }
  ")),
      tabPanel("Home", 
        fluidRow(
          column(width = 2),  # Left margin
          column(width = 8,   # Center content
            div(style = "text-align: center;",
              h2('PON-Del predictor for short protein deletions'),
              br()
            ),
            div(style = "text-align: justify; margin: 0 auto; max-width: 800px;",
              p('PON-Del is a predictor for short (1-10 amino acid) sequence retaining deletions. It was trained on an extensive set of variations and showed superior performance compared to other tools. After evaluating multiple frameworks, LightGBM was selected as the final model.')
            ),
            div(style = "text-align: center; margin-top: 30px;",
              img(src = "overview.png", height = "500px", alt = "PON-Del Overview", 
                  style = "max-width: 100%; height: auto; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1);")
            )
          ),
          column(width = 2)   # Right margin
        ),
        # Add prediction panel under home page
        fluidRow(
          column(
            width = 5,
            wellPanel(
              style = "background: #fcfcfc; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.03); padding: 20px; min-height: 500px;",
              h4('Input deletions', style = "margin-bottom: 20px;"),
              
              # Add deletion type selection
              h5('Deletion type', style = "margin-bottom: 10px; margin-top: 0;"),
              radioGroupButtons(
                inputId = "input_type",
                label = NULL,
                choices = c(
                  "Protein" = "protein",
                  "Transcript" = "transcript", 
                  "Genomic" = "genomic"
                ),
                selected = "protein",
                direction = "horizontal",
                justified = FALSE,
                size = "xs"
              ),
              
              
              textAreaInput("deletion_input", 
                label = textOutput("deletion_input_label"),
                value = protein_demo,
                placeholder = textOutput("deletion_input_placeholder"),
                rows = 10,
                width = "100%"
              ),
              div(style = "margin-top: 5px; font-size: 12px; color: #666;",
                HTML("<strong>Note:</strong> Highly recommend to use protein-level variant. If submit with transcript or genomic, we will use transvar to convert to protein level and it may miss some cases. Maximum 1000 lines allowed. Deletions starting at position 1 will return Pathogenic because deletion of the first amino acid (usually methionine) may prevent normal protein expression.")
              ),
              div(style = "margin-top: 20px;",
                actionButton("run_prediction", "Run prediction", 
                  style = "width: 100%; background-color: #4CAF50; color: white;")
              ),
              div(style = "margin-top: 20px;",
                downloadButton("download_prediction", "Download Results (CSV)", 
                  style = "width: 100%; background-color: #007bff; color: white; border: none; padding: 10px; border-radius: 4px;")
              ),
              div(style = "margin-top: 20px;",
                h5("Prediction log"),
                verbatimTextOutput("prediction_log")
              )
            )
          ),
          column(
            width = 7,
            wellPanel(
              style = "background: #fcfcfc; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.03); padding: 20px; min-height: 500px;",
              h4('Prediction results'),
              div(
                DT::dataTableOutput("prediction_table")
              ),
              # Column descriptions
              div(style = "margin-top: 15px;",
                h5("Column Descriptions:"),
                tags$ul(
                  tags$li(tags$strong("Input:"), "Original input coordinates in the format you provided"),
                  tags$li(tags$strong("RefSeq Protein:"), "Corresponding RefSeq protein identifier"),
                  tags$li(tags$strong("Deletion Start/End:"), "Protein positions of the deletion"),
                  tags$li(tags$strong("Probability:"), "Prediction score (0-1, higher = more pathogenic)"),
                  tags$li(tags$strong("Prediction:"), "P = Pathogenic (probability > 0.5), B = Benign (probability > 0.5), U = Uncertain (bootstrap P>0.05)."),
                  tags$li(tags$strong("P value:"), "Bootstrap P-value (P>0.05 = Uncertain)"),
                )
              ),
              div(style = "margin-top: 15px; font-size: 12px; color: #666; padding: 10px; background-color: #f8f9fa; border-radius: 4px;",
                HTML("<strong>Note:</strong> Deletions starting at position 1 return Pathogenic because deletion of the first amino acid (usually methionine) may prevent normal protein expression and are beyond the scope of the current prediction model. If you need binary prediction, just set probability > 0.5 as pathogenic otherwise benign, neglect the p.")
              )
            )
          )
        )
      ),
      
      tabPanel("Single Amino Acid Deletion", 
        # Title and description
        fluidRow(
          column(width = 2),  # Left margin
          column(width = 8,   # Center content
            div(style = "text-align: center; margin-bottom: 30px;",
              h2('Single amino acid deletions')
            ),
            div(style = "text-align: justify; margin: 0 auto; max-width: 800px;",
              p('This page provides precalculated results for all possible single amino acid deletions in proteins coded by MANE transcripts.'),
              p('Deletions starting at position 1 return Pathogenic because deletion of the first amino acid (usually methionine) may prevent normal protein expression.'),
              p('You can search the data in several different ways.')
            )
          ),
          column(width = 2)   # Right margin
        ),
        # First row: Selection and table
        fluidRow(
          column(
            width = 5,
            tags$style(HTML("
              .btn-group, .btn-group-justified {
                flex-wrap: wrap !important;
                display: flex !important;
              }
              .btn-group .btn, .btn-group-justified .btn {
                min-width: 160px;
                margin-bottom: 6px;
              }
            ")),
            wellPanel(
              style = "background: #fcfcfc; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.03); padding: 20px; min-height: 500px;",
              h4('Select identifier', style = "margin-bottom: 0.5em; margin-top: 20px;"),
              h5('1. Choose identifier type', style = "margin-bottom: 0.5em; margin-top: 20px;"),
              radioGroupButtons(
                inputId = "id_type",
                label = NULL,
                choices = id_types,
                selected = "NP_id",
                direction = "horizontal",
                justified = FALSE,
                size = "xs"
              ),
              h5('2. Enter ID', style = "margin-bottom: 0.5em; margin-top: 20px;"),
              selectizeInput("protein_id", label = NULL,
                choices = NULL,
                options = list(placeholder = "Type to search..."),
                width = "100%"),
              div(style = "margin-top: 25px;",
                downloadButton("download_plot", "Download Heat Map (PNG)", 
                  style = "width: 100%; margin-bottom: 10px; background-color: #28a745; color: white; border: none; padding: 8px; border-radius: 4px;"),
                downloadButton("download_table", "Download Results (CSV)", 
                  style = "width: 100%; background-color: #007bff; color: white; border: none; padding: 8px; border-radius: 4px;")
              ),
              tags$div(
                style = "font-size: 13px; color: #888; margin-top: 6px; margin-bottom: 10px;",
                HTML('<br> Note:<br> PON-Del is developed based on MANE selectedRefSeq protein identifiers. <br><br> If you select a different identifier type, the corresponding RefSeq protein will be displayed. <br><br> The <strong>one-to-one mapping</strong> is defined by the <a href="https://ftp.ncbi.nih.gov/refseq/MANE/MANE_human/release_1.4/" target="_blank">MANE v1.4</a>.')
              )
            )
          ),
          column(
            width = 7,
            wellPanel(
              style = "background: #fcfcfc; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.03); padding: 20px; min-height: 500px;",
              h4('Predicted deletion pathogenicity'),
              div(
                DT::dataTableOutput("protein_table")
              ),
              tags$div(
                style = "font-size: 13px; color: #888; margin-top: 6px; margin-bottom: 10px;",
                htmlOutput("prediction_summary")
              )
            )
          )
        ),
        # Second row: Plot (reduce top margin)
        fluidRow(
          column(width = 12,
            wellPanel(
              style = "background: #fcfcfc; borRr-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.03); padding: 20px; margin-top: -10px;",
              h4('Heatmap for predicted pathogenicity', style = "margin-bottom: 10px;"),
              imageOutput("protein_plot", height = "auto")
            )
          )
        )
      ),

      tabPanel("About", 
        fluidRow(
          column(width = 2),  # Left margin
          column(width = 8,   # Center content
            div(
              style = "max-width: 800px; margin: 0 auto; text-align: left;",
              h2('About', style = "text-align: center;"),
              h3('Data and Code', style = "margin-top: 30px;"),
              p('The datasets used for training and testing the tool are available ',
                a('here', href = 'data_pondel.csv', target = '_blank'), '.'),
              p('The training data and code are available on ', 
                a('GitHub', href = 'https://github.com/zhanghaoyang0/pon_del_public', target = '_blank'), '.'),
              
              h3('Citing PON-Del', style = "margin-top: 30px;"),
              p('A manuscript describing the predictor has been submitted. In the meantime use URL for citation.'),
              
              h3('Contact', style = "margin-top: 30px;"),
              p(
                'The tool was developed by Haoyang Zhang and Muhammad Kabir, supervised by Mauno Vihinen.'
              ),
              p(  
                'If you have any problems, please contact Haoyang (haoyang.zhang@med.lu.se) or Mauno (mauno.vihinen@med.lu.se).'
              )
            )
          ),
          column(width = 2)   # Right margin
        )
      )
    )
)

server = shinyServer(function(input, output, session) {
  # Update protein ID choices based on selected ID type
  observe({
    selected_type <- input$id_type
    choices <- unique(id_mapping[[selected_type]])
    choices <- choices[!is.na(choices)]
    updateSelectizeInput(session, "protein_id",
                        choices = choices,
                        selected = NULL,
                        server = TRUE)
  })

  # Dynamic labels and placeholders for deletion input based on type
  output$deletion_input_label <- renderText({
    input_type <- input$input_type
    switch(input_type,
           "protein" = "Enter deletions. Format: NP_id,deletion start(protein),deletion end(protein)",
           "transcript" = "Enter deletions. Format: ENST_id,deletion start(transcript),deletion end(transcript)",
           "genomic" = "Enter deletions. Format: chr,deletion start(genomic),deletion end(genomic)",
           "Enter deletions. Format: NP_id,deletion start(protein),deletion end(protein)")
  })
  
  output$deletion_input_placeholder <- renderText({
    input_type <- input$input_type
    switch(input_type,
           "protein" = "NP_000007.1,115,116\nNP_000008.1,104,104\nNP_000009.1,299,299",
           "transcript" = "NM_000007.3,345,348\nNM_000008.2,312,312\nNM_000009.3,897,897",
           "genomic" = "chr1:69091-69093\nchr2:476302-476304\nchr3:123456-123458",
           "NP_000007.1,115,116\nNP_000008.1,104,104\nNP_000009.1,299,299")
  })

  # Get NP_id from selected ID
  get_np_id <- reactive({
    req(input$protein_id, input$id_type)
    selected_type <- input$id_type
    selected_id <- input$protein_id
    
    if (selected_type == "NP_id") {
      return(selected_id)
    }
    
    np_id <- id_mapping$NP_id[id_mapping[[selected_type]] == selected_id]
    if (length(np_id) == 0) return(NULL)
    return(np_id[1])
  })

  # Update text area with demo input when deletion type changes
  observeEvent(input$input_type, {
    demo_text <- switch(input$input_type,
                        "protein" = protein_demo,
                        "transcript" = transcript_demo,
                        "genomic" = genomic_demo,
                        "")
    updateTextAreaInput(session, "deletion_input", value = demo_text)
  })

  output$protein_plot <- renderImage({
    req(get_np_id())
    np_id <- get_np_id()
    plot_path <- file.path("/srv/shiny-server/db/pon_del/plot/1aa_del", paste0(np_id, ".png"))
    
    if (!file.exists(plot_path)) {
      return(list(
        src = "",
        contentType = "image/png",
        width = "100%",
        height = "auto"
      ))
    }
    
    list(
      src = plot_path,
      contentType = "image/png",
      width = "100%",
      height = "auto"
    )
  }, deleteFile = FALSE)

  # Download handler for plot
  output$download_plot <- downloadHandler(
    filename = function() {
      np_id <- get_np_id()
      if (is.null(np_id) || np_id == "") {
        return("heatmap.png")
      }
      paste0(np_id, "_heatmap.png")
    },
    content = function(file) {
      np_id <- get_np_id()
      if (is.null(np_id) || np_id == "") {
        # Return empty file if no protein selected
        writeLines("No protein selected", file)
        return()
      }
      
      plot_path <- file.path("/srv/shiny-server/db/pon_del/plot/1aa_del", paste0(np_id, ".png"))
      if (file.exists(plot_path)) {
        file.copy(plot_path, file)
      } else {
        # Return error message if file doesn't exist
        writeLines("Heat map file not found", file)
      }
    }
  )

  # Summary statistics for predictions
  output$prediction_summary <- renderText({
    req(get_np_id())
    np_id <- get_np_id()
    table_path <- file.path("/srv/shiny-server/db/pon_del/data/1aa_del/output", paste0(np_id, ".csv"))
    
    if (!file.exists(table_path)) {
      return("")
    }
    
    df <- read.csv(table_path, check.names = TRUE, strip.white = TRUE)
    if (nrow(df) == 0) {
      return("")
    }
    
    # Map labels: 1 = Pathogenic, 0 = Benign, 2 = Uncertain
    n_P <- sum(df$pondel_label == 1, na.rm = TRUE)
    n_B <- sum(df$pondel_label == 0, na.rm = TRUE)
    n_U <- sum(df$pondel_label == 2, na.rm = TRUE)
    n_total <- nrow(df)
    
    if (n_total > 0) {
      HTML(paste0(
        "<br> Summary:<br> Total: ", n_total, 
        " | Pathogenic (P): ", n_P, 
        " | Benign (B): ", n_B, 
        " | Uncertain (U): ", n_U
      ))
    } else {
      ""
    }
  })
  
  output$protein_table <- DT::renderDataTable({
    req(get_np_id())
    np_id <- get_np_id()
    table_path <- file.path("/srv/shiny-server/db/pon_del/data/1aa_del/output", paste0(np_id, ".csv"))
    
    if (!file.exists(table_path)) {
      return(data.frame(Message = "No data available for this protein."))
    }
    
    df <- read.csv(table_path, check.names = TRUE, strip.white = TRUE)
    if (nrow(df) == 0) {
      data.frame(Message = "No data rows in CSV.")
    } else {
      # Map labels: 1 = Pathogenic, 0 = Benign, 2 = Uncertain
      df$pondel_label <- ifelse(df$pondel_label == 1, "P",
                                ifelse(df$pondel_label == 0, "B",
                                       ifelse(df$pondel_label == 2, "U", "N/A")))
      
      # Rename columns for display
      df = df%>%dplyr::rename(
        "RefSeq Protein" = NP_id,
        "Probability" = pondel_pred,
        "Prediction" = pondel_label,
        "Deletion Position" = start
      )
      # Round probability and map label
      df$`Probability` <- round(df$`Probability`, 2)

      DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
    }
  })

  # Download handler for table
  output$download_table <- downloadHandler(
    filename = function() {
      np_id <- get_np_id()
      if (is.null(np_id) || np_id == "") {
        return("predictions.csv")
      }
      paste0(np_id, "_predictions.csv")
    },
    content = function(file) {
      np_id <- get_np_id()
      if (is.null(np_id) || np_id == "") {
        # Return empty file if no protein selected
        write.csv(data.frame(Message = "No protein selected"), file, row.names = FALSE)
        return()
      }
      
      table_path <- file.path("/srv/shiny-server/db/pon_del/data/1aa_del/output", paste0(np_id, ".csv"))
      if (file.exists(table_path)) {
        df <- read.csv(table_path, check.names = TRUE, strip.white = TRUE)
        if (nrow(df) > 0) {
          # Rename columns for download
          df = df%>%dplyr::rename(
            "RefSeq Protein" = NP_id,
            "Probability" = pondel_pred,
            "Prediction" = pondel_label,
            "Deletion Position" = start
          )
          # Round probability and map label
          df$`Probability` <- round(df$`Probability`, 2)
          write.csv(df, file, row.names = FALSE)
        } else {
          write.csv(data.frame(Message = "No data available"), file, row.names = FALSE)
        }
      } else {
        write.csv(data.frame(Message = "File not found"), file, row.names = FALSE)
      }
    }
  )

  # Handle prediction tab
  prediction_results <- reactiveVal(NULL)
  prediction_log <- reactiveVal("")
  is_running <- reactiveVal(FALSE)
  
  # Simple logging function that directly updates UI
  log_message <- function(msg) {
    # Get current time
    time <- format(Sys.time(), "%H:%M:%S")
    # Format message with time
    formatted_msg <- sprintf("[%s] %s\n", time, msg)
    # Update the reactive value
    prediction_log(paste0(prediction_log(), formatted_msg))
    # Force UI update
    invalidateLater(100)
  }
  
  # Update button text and state based on running status
  observe({
    if (is_running()) {
      updateActionButton(session, "run_prediction", 
                        label = "Running...", 
                        icon = icon("spinner", class = "fa-spin"))
    } else {
      updateActionButton(session, "run_prediction", 
                        label = "Run Prediction", 
                        icon = NULL)
    }
  })
  
  observeEvent(input$run_prediction, {
    # Set running state
    is_running(TRUE)
    
    # Clear previous log and results
    prediction_log("")
    prediction_results(NULL)
    
    # Start logging
    log_message("Starting prediction...")
    
    # Log the selected deletion type
    input_type <- input$input_type
    log_message(sprintf("Deletion type: %s", input_type))
    
    # Validate input
    if (is.null(input$deletion_input) || input$deletion_input == "") {
      log_message("Error: No input data provided")
      is_running(FALSE)
      return()
    }
    
    # Check line limit (maximum 1000 lines)
    input_lines <- strsplit(input$deletion_input, "\n")[[1]]
    if (length(input_lines) > 1000) {
      log_message(sprintf("Error: Too many lines (%d). Maximum 1000 lines allowed", length(input_lines)))
      is_running(FALSE)
      return()
    }
    
    log_message(sprintf("Processing %d lines of input data", length(input_lines)))
    
    # Create temp directory with full path
    temp_dir <- file.path(getwd(), "temp")
    if (!dir.exists(temp_dir)) {
      dir.create(temp_dir, recursive = TRUE, mode = "0777")
    }

    if (file.access(temp_dir, 2) != 0) {
      log_message("Warning: temp directory is not writable")
    }

    # Generate task ID
    task_id <- format(Sys.time(), "%Y%m%d%H%M%S")
    input_file <- file.path(temp_dir, paste0(task_id, ".csv"))
    
    # Write input data with error handling
    tryCatch({
      writeLines(input$deletion_input, input_file)
      if (!file.exists(input_file)) {
        log_message("Error: Input file was not created")
        is_running(FALSE)
        return()
      }
    }, error = function(e) {
      log_message(sprintf("Error saving input: %s", e$message))
      is_running(FALSE)
      return()
    })

    # Run pipeline with deletion type parameter
    log_message("Running prediction pipeline...")
    pipeline_cmd <- paste("sh code/pipeline.sh", task_id, input$input_type)

    tryCatch({
      # Run pipeline and capture output
      result <- system(pipeline_cmd, intern = TRUE)
      log_message("Pipeline completed")
    }, error = function(e) {
      log_message(sprintf("Pipeline error: %s", e$message))
      is_running(FALSE)
      return()
    })

    # Check results
    results_file <- file.path(temp_dir, paste0(task_id, "_predictions.csv"))

    if (file.exists(results_file)) {
      tryCatch({
        results <- read.csv(results_file)
        prediction_results(results)
        log_message(sprintf("Prediction completed successfully! Found %d results", nrow(results)))
      }, error = function(e) {
        log_message(sprintf("Error reading results: %s", e$message))
      })
    } 

    log_message("Prediction finished")
    
    # Reset running state
    is_running(FALSE)
  })
  
  # Update the UI to show logs
  output$prediction_log <- renderText({
    req(prediction_log())
    prediction_log()
  })
  
  # Helper function to format prediction data
  format_prediction_data <- function(df, for_display = TRUE) {
    # Handle special cases and format the results
    # Position 1 deletions are set to "P" (Pathogenic) by the Python script
    # Handle both numeric and string values for pondel_pred
    df$pondel_pred <- ifelse(df$pondel_pred == "P", "P", 
                            ifelse(is.na(df$pondel_pred) | is.nan(df$pondel_pred), "N/A", 
                                   round(as.numeric(df$pondel_pred), ifelse(for_display, 2, 3))))
    
    # Helper function to format labels consistently
    format_label <- function(x) {
      ifelse(x == "P", "P", 
             ifelse(x == "B", "B", 
                    ifelse(x == "U", "U", "N/A")))
    }
    
    # Apply label formatting to all label columns
    df$pondel_label <- format_label(df$pondel_label)
    
    # Helper function to format numeric columns
    format_numeric <- function(x, for_display = TRUE) {
      ifelse(x == "N/A", "N/A", 
             ifelse(is.na(x) | is.nan(x), "N/A", 
                    round(as.numeric(x), ifelse(for_display, 2, 4))))
    }
    
    # Apply numeric formatting for p-value
    if ("pondel_pred_pvalue" %in% colnames(df)) {
      df$pondel_pred_pvalue <- format_numeric(df$pondel_pred_pvalue, for_display)
    }
    
    return(df)
  }
  
  output$prediction_table <- DT::renderDataTable({
    req(prediction_results())
    df <- prediction_results()
    
    # Format data for display
    df <- format_prediction_data(df, for_display = TRUE)
    
    # Remove pondel_pred_label column if it exists
    if ("pondel_pred_label" %in% colnames(df)) {
      df <- df %>% dplyr::select(-pondel_pred_label)
    }
    
    # Move Input column to first position if it exists
    if ("Input" %in% colnames(df)) {
      df <- df %>% dplyr::select(Input, everything())
      print("Moved Input column to first position in prediction table")
    }
    
    # Rename columns for display
    df <- df %>% dplyr::rename(
      "RefSeq Protein" = NP_id,
      "Probability" = pondel_pred,
      "Prediction" = pondel_label,
      "Deletion Start" = start,
      "Deletion End" = end
    )
    
    # Add optional columns if they exist
    if ("Input" %in% colnames(df)) {
      df <- df %>% dplyr::rename("Input" = Input)
    }
    if ("pondel_pred_pvalue" %in% colnames(df)) {
      df <- df %>% dplyr::rename("P value" = pondel_pred_pvalue)
    }
    
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$download_prediction <- downloadHandler(
    filename = function() {
      paste0("pondel_prediction_results_", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv")
    },
    content = function(file) {
      req(prediction_results())
      df <- prediction_results()
      
      # Format data for download (3 decimal places)
      download_df <- format_prediction_data(df, for_display = FALSE)
      
      # Remove pondel_pred_label column if it exists
      if ("pondel_pred_label" %in% colnames(download_df)) {
        download_df <- download_df %>% dplyr::select(-pondel_pred_label)
      }
      
      # Move Input column to first position if it exists
      if ("Input" %in% colnames(download_df)) {
        download_df <- download_df %>% dplyr::select(Input, everything())
        print("Moved Input column to first position in download")
      }
      
      # Rename columns for download (with underscores)
      if ("Input" %in% colnames(download_df)) {
        colnames(download_df)[colnames(download_df) == "Input"] <- "Input_Coordinates"
      }
      colnames(download_df)[colnames(download_df) == "NP_id"] <- "RefSeq_Protein"
      colnames(download_df)[colnames(download_df) == "pondel_pred"] <- "Probability"
      colnames(download_df)[colnames(download_df) == "pondel_label"] <- "Prediction"
      if ("pondel_pred_pvalue" %in% colnames(download_df)) {
        colnames(download_df)[colnames(download_df) == "pondel_pred_pvalue"] <- "P_value"
      }
      colnames(download_df)[colnames(download_df) == "start"] <- "Deletion_Start"
      colnames(download_df)[colnames(download_df) == "end"] <- "Deletion_End"
      
      write.csv(download_df, file, row.names = FALSE)
    }
  )
})

shinyApp(ui, server) 