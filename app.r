# ====================================================================
# prepare
# ====================================================================
library(shiny)
library(shinydashboard)
library(DT)
library(digest)


ui = shinyUI(navbarPage(
  title = "PON-mt-tRNA",
      tabPanel("Home", 
        mainPanel(
          h3('PON-mt-tRNA'),
          'PON-mt-tRNA is a posterior probability-based tool for classification of human mitochondrial tRNA (mt-tRNA) variations. ', 
          br(), br(), 'It predicts the probability of pathogenicity using 20 machine learning (ML) predictors. The predicted probability is used as 
          prior probability to integrate the evidence submitted by the user. ', 
          br(), br(), 'If no evidence is submitted, PON-mt-tRNA classifies the variations based on ML predicted probability. 
          The ML predictors utilize features representing evolutionary conservation, sequence context, secondary structure and tertiary interactions.',  
          br(), br(), 'Balanced accuracy and Matthews correlation coefficient (MCC) for integrated classifier based on posterior probability are 0.99 and 0.95, 
          respectively and for ML predictor based probability of pathogenicity are 0.81 and 0.56, respectively on independent test dataset.',  
          br(), br(), 'PON-mt-tRNA is trained using variations that were classified previously by using evidence-based classification method by Yarham et al.', 
          a(href = "http://www.ncbi.nlm.nih.gov/pubmed/21882289", 'Paper link'),
          h3('Citation'),
          a(href = "https://www.gu.se/en/about/find-staff/abhishekniroula", 'Abhishek Niroula'),
          ' and ', 
          a(href = "https://portal.research.lu.se/en/persons/mauno-vihinen", 'Mauno Vihinen.'),
          'PON-mt-tRNA: a multifactorial probability-based method for classification of mitochondrial tRNA variations.'
          , em('Nucleic Acids Res'), '. (2016) 44(5): 2020-2027.', 
          a(href = "http://nar.oxfordjournals.org/content/early/2016/02/02/nar.gkw046.abstract", 'Paper link'),
          h3('Comment and feedback'),
          'The server is maintain by Haoyang Zhang. If you meet have questions, please contact haoyang.zhang@med.lu.se',
          br(), br(), img(src = "pon_mtRNA_logo.png",   height = "150px", width = "150px", align = "center",  alt = "logo"))
      ),
      tabPanel("Datasets", 
        mainPanel(
          h3('Datasets used in PON-mt-tRNA'),
          'The datasets described in the paper are available for download.', 
          'It predicts the probability of pathogenicity using 20 machine learning (ML) predictors. The predicted probability is used as 
          prior probability to integrate the evidence submitted by the user. ', 
          h3('Feature matrix for training and test data'),
          'This file contains the feature matrix used for developing PON-mt-tRNA. It contains 91 pathogenic and 55 neutral variations. 40 pathogenic 
          and 40 neutral variations were selected by random sampling without replacement for training and the remaining variations were used for testing the method.',  
          a(href = "PON-mt-tRNA_feature_matrix.xlsx", 'Download'),
          h3('Additional variation dataset'),
          'This file contains predictions of PON-mt-tRNA for additional variants obtained from MITOMAP, mtDB and mtSNP databases. 
          The variations that were present in the PON-mt-tRNA training and test dataset were excluded from this dataset.', 
          a(href = "PON-mt-tRNA_additional_datasets.xlsx", 'Download'), 
          h3('PON-mt-tRNA predictions'),
          'This file contains predictions of PON-mt-tRNA for all possible single nucleotide substitutions at each position in the 22 human mt-tRNA. 
          The classification is based on ML predicted probability of pathogenicity.', 
          a(href = "PON-mt-tRNA_predictions.txt", 'Download')
        )
      ), 
      tabPanel("Input and output", 
        mainPanel(
          h3('Submit queries to PON-mt-tRNA'), 
          'PON-mt-tRNA requires variation location, and the reference and altered nucleotides separated by comma (csv). In addition, the users can submit evidence of segregation, 
          biochemical test and histochemical test which are optional. The evidence should be submitted in the following manner:',
          h4('Example input'), 
          "POS,REF,ALT,Segregation,Biochemistry,Histochemistry",
          br(), "7505,A,G",br(), "7505,A,G,1", br(), "8300,T,C,NA,1,NA",br(), "7316,T,C,0,NA,1", 
          br(), "7316,T,A,0,NA,0",br(), "7316,T,C,1,1,1",br(), "7316,T,A,NA,NA,NA",br(), "7316,T,C,0,0,0",
          h3('Note'),
          h4('POS, REF, ALT'),
          'mtDNA location, Reference nucleotide, Altered nucleotide',
          h4('Segregation'),
          '1: Segregation of variation with disease',
          br(), '0: No segregation with disease',
          br(), 'NA: Not available',
          h4('Biochemistry'),
          '1: Biochemical defect in complexes I, III or IV',
          br(), '0: No biochemical defect in the complexes I, III and IV',
          br(), 'NA: Not available',
          h4('Histochemistry'),
          '1: Histochemical evidence of mitochondrial disease',
          br(), '0: No histochemical evidence of mitochondrial disease',
          br(), 'NA: Not available',
          br(), br(), 'The evidence should be submitted in the order Segregation, Biochemistry and Histochemistry.',
          br(), 'If evidence of Segregation is not known but others are known, NA should be used for segregation and then other should be provided.', 
          br(), 'Multiple variations can be submitted in a single query. PON-mt-tRNA only classifies the single nucleotide substitutions and therefore 
          does not accept other types of variations as input. Alternatively, a file containing the variations in the same format as described above can be uploaded.',
          h3('PON-mt-tRNA output'),
          'The result is like:',
          br(), 'mt-tRNA	Variation	ML_probability_of_pathogenicity	Evidence	Posterior_probability_of_pathogenicity	Classification',
          br(), 'Ser(UCN)	7505,A,G	0.46	NA, NA, NA		Likely neutral',
          br(), 'Ser(UCN)	7505,A,G	0.46	1, NA, NA	0.66	VUS',
          br(), 'Lys	8300,T,C	0.45	NA, 1, NA	0.87	VUS',
          br(), 'Ala	5644,T,C	0.56	0, NA, 1	0.56	VUS',
          br(), 'Ala	5644,T,A	0.56	0, NA, 0	0.01	Likely neutral',
          br(), 'Ala	5644,T,C	0.56	1, 1, 1	0.99	Pathogenic',
          br(), 'Pro	15984,A,T	0.52	NA, NA, NA		Likely pathogenic',
          br(), 'Pro	15984,A,C	0.52	0, 0, 0	0.0	Neutral',
          br(), br(), 'The result file contains the following contents:',
          br(), '1. mt-tRNA affected by the variation', 
          br(), '2. Variation',
          br(), '3. ML probability of pathogenicity based on 20 ML predictors. It ranges from 0 to 1.',
          br(), '4. Evidence submitted by the user.',
          br(), '5. Posterior probability of pathogenicity after combining the ML predicted probability of pathogenicity and evidence submitted by the user. It ranges from 0 to 1.',
          br(), '6. Classification of variation. The classification is based on posterior probability if there is a posterior probability for the variation in column 5. Otherwise, the classification is based on ML predicted probability of pathogenicity in column 3.'
        )
      ), 

      tabPanel("Submit queries", 
        fluidRow(
          column(width = 5,
            mainPanel(h3('Submit queries'),
              h4('Submit with text'),
              tags$textarea(id = "text_input", rows = 15, style = "width:100%; resize:both;", 
                "POS,REF,ALT,Segregation,Biochemistry,Histochemistry\n7505,A,G\n7505,A,G,1\n8300,T,C,NA,1,NA\n5644,T,C,0,NA,1\n5644,T,A,0,NA,0\n5644,T,C,1,1,1\n15984,A,T,NA,NA,NA\n15984,A,C,0,0,0"),
                h4('Or upload a file containing variations:'),
              fileInput("table_input", ""),
              a(href = "example.csv", 'Example'),
          )),
          column(width = 7,
          mainPanel(h3('Prediction result'),
            DT::dataTableOutput("table_pred"),
            downloadButton("download_pred", "Download")
            )
          )
        )
      ),

      tabPanel("Disclaimer", 
        mainPanel(
          h3('Disclaimer'),
          'This non-profit server, its associated data and services are for research purposes only. The responsibility of Protein Structure 
          and Bioinformatics Group, Lund University is limited to applying the best efforts in providing and publishing good programs and data. 
          The developers have no responsibility for the usage of results, data or information which this server has provided.',
          h3('Liability'),
          'In preparation of this site and service, every effort has been made to offer the most current and correct information possible.
          We render no warranty, express or implied, as to its accuracy or that the information is fit for a particular purpose, and will not 
          be held responsible for any direct, indirect, putative, special, or consequential damages arising out of any inaccuracies or 
          omissions. In utilizing this service, individuals, organizations, and companies absolve Lund University or any of their employees 
          or agents from responsibility for the effect of any process, method or product or that may be produced or adopted by any part, 
          notwithstanding that the formulation of such process,
          method or product may be based upon information provided here.'
        )
      )
    )
  )



server = shinyServer(function(input, output) {
  data <- reactive({
    jobid <- digest(paste0(Sys.time(), sample(1:9999, 1)), algo = "md5")
    file_query <- sprintf('temp/%s.query', jobid)
    file_proc <- sprintf('temp/%s_proc_evid.txt', jobid) # temp file
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
    file.remove(file_query, file_proc, file_pred)
    pred
  })
  
  output$table_pred <- DT::renderDataTable({
    data()
  }, options = list(pageLength = 10, scrollX = TRUE))
  

})


shinyApp(ui, server)