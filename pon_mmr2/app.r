# ====================================================================
# prepare
# ====================================================================
library(shiny)
library(shinydashboard)
library(DT)
library(digest)


ui = shinyUI(fluidPage(
    tags$head(tags$title("PON-MMR2")),
      mainPanel(
        br(), br(), img(src = "pon_mmr_logo.png", align = "center",  alt = "logo"),
        h3('About PON-MMR2'),
        'PON-MMR2 is a machine-learning based tool that classifies amino acid substitutions in mismatch repair (MMR) proteins: MLH1, MSH2, MSH6 and PMS2.
        It utilizes random forest algorithm and features from amino acid features and evolutionary information. ', 
        br(), br(), 'It is trained and tested on variants obtained from', 
        a(href = "https://www.insight-group.org/variants/databases/", 'InSiGHT database'),
        'that are classified into benign (Classes 1 and 2) and harmful (Classes 4 and 5) (',  
        a(href = "https://pubmed.ncbi.nlm.nih.gov/24362816/", 'Thompson et al.'), ')',
        ', and, variants obtained from',  
        a(href = "https://structure-next.med.lu.se/VariBench/index.php", 'VariBench database'),
        'that were used to train PON-MMR (',
        a(href = "https://pubmed.ncbi.nlm.nih.gov/22290698/", 'Ali et al.'), ').',
        'Balanced accuracy and Matthews correlation coefficient (MCC) of PON-MMR2 are 0.89 and 0.78, respectively in leave-one-out cross-validation and 0.85 and 0.67, respectively on an independent test dataset.',
        h3('PON-MMR2 predictions'), 
        'The reference protein sequences used in PON-MMR2 are obtained from UniProtKB. You can find them in the links below.', 
        br(), a(href = "MLH1.fa", 'MLH1|P40692'), 
        br(), a(href = "MSH2.fa", 'MSH2|P43246'), 
        br(), a(href = "MSH6.fa", 'MSH6|P52701 '), 
        br(), a(href = "PMS2.fa", 'PMS2|P54278'), 
        br(), br(),'We predicted all possible amino acid substitutions at each position in the four MMR proteins using PON-MMR2. You can find them in the link below.',
        br(), a(href = "PON-MMR2_predictions.txt", 'Download all predictions'), 
        br(), 'The result file contains the following contents:',
        br(), '1. Gene name(s)',
        br(), '2. UniProtKB accession number(s) for the reference protein sequence',
        br(), '3. Amino acid substitution(s)',
        br(), '4. Positions of amino acid substitution in the reference protein sequence',
        br(), '5. Original amino acid in the reference protein sequence',
        br(), '6. New amino acid that substitutes the original amino acid',
        br(), '7. Probability of pathogenicity predicted by PON-MMR2. It ranges from 0 (benign) to 1 (harmful).',
        br(), '8. Classification of the variation as pathogenic or neutral based on the probability of pathogenicity (Column 7).',
        br(), '9. InSiGHT Class. It is provided for the amino acid substitutions that are present in PON-MMR2 training and test datasets.', br(), br(),
        h4('Query PON-MMR2 predictions'), 
        DT::dataTableOutput("table_pred"),
        h3('Reference'),
        a(href = "https://www.gu.se/en/about/find-staff/abhishekniroula", 'Abhishek Niroula'),
        ' and ', 
        a(href = "https://portal.research.lu.se/en/persons/mauno-vihinen", 'Mauno Vihinen.'),
        br(), 'Classification of amino acid substitutions in mismatch repair proteins using PON-MMR2.'
        , em('Hum Mutat'), '.2015.', 
        a(href = "https://onlinelibrary.wiley.com/doi/full/10.1002/humu.22900", 'Paper link'),
        h3('Comment and feedback'),
        'The server is maintain by Haoyang Zhang. If you meet have questions, add an issue at',
        a(href = "https://github.com/zhanghaoyang0/lu_psb_rshiny_server", 'our Github'),
        'or send an email to haoyang.zhang@med.lu.se',
      )
    )
  )
  




server = shinyServer(function(input, output) {
  data <- reactive({
    jobid <- digest(paste0(Sys.time(), sample(1:9999, 1)), algo = "md5")
    file_query <- sprintf('temp/%s.query', jobid)
    file_pred <- sprintf('temp/%s_results.txt', jobid)
    if (!is.null(input$table_input)) {
      df <- read.csv(input$table_input$datapath, header = TRUE)
    } else {
      df <- read.csv(text = input$text_input, header = TRUE)
    }
    write.csv(df, file_query, row.names = FALSE, quote=F)
    command = sprintf('python3 query.py %s', jobid)
    system(command)
    pred = read.table(file_pred, sep='\t', header=T)
    file.remove(file_query, file_pred)
    pred
  })
  
  output$table_pred <- DT::renderDataTable({
    read.delim('www/PON-MMR2_predictions.txt')
  }, options = list(pageLength = 10, scrollX = TRUE))
  

})


shinyApp(ui, server)