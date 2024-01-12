# ====================================================================
# prepare
# ====================================================================
library(shiny)
library(shinydashboard)
library(DT)
library(digest)


ui <- fluidPage(
    tags$head(tags$title("PON-BT")),
    h1('PON-BTK'),
    'PON-BTK is a method for classifying variations in the kinase domain of Bruton tyrosine kinase (BTK) related to X-linked agammaglobulinemia (XLA) to disease-causing and harmful. 
    The method is using the semi-supervised classification method', 
    a(href = "https://cran.r-project.org/web/packages/upclass/index.html", "'upclass'"), 
    '. It uses an expectation maximization (EM) algorithm to obtain maximum likelihood estimates of the model parameters and classifications for the unlabeled data.',
    br(), br(), 
    'The PON-P probabilities of pathogenicity were obtained directly from PON-P output and the consensus pathogenicity is the proportion of pathogenic predictions of the methods having Matthews correlation coefficient (MCC) > 0.5. The precalculated results (see below) include in addition to the PON-BTK prediction also other information about the variants.', 
    h3('Query PON-BTK predictions'),
    DT::dataTableOutput("table_pred"),
    downloadButton("download_pred", "Download all PON-BTK predictions"),
    h3('Reference'),
    'Jouni Väliaho, Imrul Faisal, Csaba Ortutay, C. I. Edvard Smith and',
    a(href = "https://portal.research.lu.se/en/persons/mauno-vihinen", 'Mauno Vihinen.'),
    br(), 'Characterization of all possible single nucleotide change-caused amino acid substitutions in the kinase domain of Bruton tyrosine kinase.'
    , em(' Hum Mutat'), '. 2015.', 
    a(href = "https://onlinelibrary.wiley.com/doi/full/10.1002/humu.22791", 'Paper link'),
    h3('Comment and feedback'),
    'The server is maintain by Haoyang Zhang. If you have questions, add an issue at',
    a(href = "https://github.com/zhanghaoyang0/pon_btk", 'Github'),
    'or send an email to haoyang.zhang@med.lu.se',
    br(), br(), img(src = "pon_btk_logo.png", align = "center",  alt = "logo")
  )


server <- function(input, output) {
  output$table_pred <- DT::renderDataTable({
    read.delim('www/pon_btk_predictions.csv')
  })

}

shinyApp(ui, server)