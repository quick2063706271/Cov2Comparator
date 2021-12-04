library(shiny)
library(shinybusy)

pdf(NULL)

ui <- fluidPage(
  tags$head(tags$style(
    HTML('
         #sidebar {
            background-color: #F8F8FF;
        }

        body, label, input, button, select {
          font-family: "Arial";
        }')
  )),
  titlePanel(tags$p("Cov2Comparator")),
  helpText("Select your countries or upload fasta files to compare genome"),
  sidebarLayout(
    sidebarPanel(
      id = "sidebar",
      tags$p(
        "Cov2Comparator can analyze SARS covid-2 genome across different
        geographic regions including performing multiple sequence alignments
        and phylogenetic tree."
      ),
      tags$p(
        "Select countries of your interests or upload your own fasta files to
        compare"
      ),
      selectInput(inputId = "selectRegion",
                  label = "Select regions of your interest",
                  choices = c("Wuhan" = "Wuhan",
                              "Canada" = "Canada",
                              "Italy" = "Italy",
                              "USA" = "USA",
                              "Kenya" = "Kenya",
                              "Bahrain" = "Bahrain",
                              "Germany" = "Germany",
                              "Pakistan" = "Pakistan",
                              "Britain" = "Britain",
                              "Thailand" = "Thailand"
                  ),
                  selected = c("Wuhan", "Canada"),
                  multiple = TRUE
      ),
      tags$p(
        "Upload your own genome sequences."
      ),
      # Select country/genome (selectInput)
      fileInput(inputId = "fastaFile", label = "Upload your Fasta file",
                accept = c(
                  ".fasta"
                )
      ),
      fileInput(inputId = "nameToRegionFile", label = "Upload nameToRegion file",
                accept = c(
                  ".txt",
                  ".csv"
                )
      ),
      actionButton(inputId = "upload", label = "Upload Genome"),
      br(),
      br(),
      selectInput(inputId = "selectAlgorithm",
                  label = "Select algorithm to run multiple sequence alignment",
                  choices = c("ClustalW" = "clustalw",
                              "Muscle" = "muscle"
                  ),
                  selected = c("ClustalW"),
                  multiple = FALSE
      ),
      actionButton(inputId = "msa", label = "Run Comparison!"),
      hr(),
      titlePanel("Modify Alignment Plot"),
      tags$p("Select the range for alignment: "),
      fluidRow(column(5,
                      textInput(inputId = "startIdx",
                                label = "From:",
                                value = 1,
                                width = "100px")
      ),
      column(5, ofset = 3,
             textInput(inputId = "endIdx",
                       label = "To:",
                       value = 1,
                       width = "100px")
      )),
      hr(),
      titlePanel("Modify Tree Plot"),
      tags$p("Select if you want to show region in the plot"),
      checkboxInput("showRegion",
                    label = "show region",
                    value = TRUE),
      selectInput(inputId = "selectTreeType",
                  label = "Select type of tree you want to plot",
                  choices = c("phylogram" = "phylogram",
                              "cladogram" = "cladogram",
                              "fan" = "fan",
                              "unrooted" = "unrooted",
                              "radial" = "radial"
                  ),
                  selected = c("phylogram"),
                  multiple = FALSE
      ),
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("MSA heatMap",
                           plotOutput("plotMsa"),
                           shinybusy::add_busy_spinner(spin = "self-building-square",
                                                       color = "#51A0E9",
                                                       position = "bottom-right"),
                           width = "auto",
                           height = "600px"
                           ),
                  tabPanel("Plot Phylo Tree",
                           plotOutput("plotTree")
                           )
      )
    )
  )
)

server <- function(input, output, session) {
  stringSet <- reactiveValues(data = NULL)
  # Upload files
  observeEvent(input$upload,{
    if (is.null(input$fastaFile)) {
      updateActionButton(session, "upload", label = "Upload Failed! No files uploaded")
    } else {
      uploadedSet <- readGenome(fastaFile = input$fastaFile$datapath,
                                nameToRegionsFile = input$nameToRegionFile$datapath)
      stringSet$data <- uploadedSet
      updateActionButton(session, "upload", label = "Upload Success!!!")
    }
  })
  alignment <- reactiveValues(data = NULL)
  plotTwoGraphs <- reactive({
    output$plotMsa <- renderPlot({
        plotAlignment(alignment$data,
                      startIdx = as.numeric(input$startIdx),
                      endIdx = as.numeric(input$endIdx))
    })
    output$plotTree <- renderPlot({
      plotTree(tree = alignment$tree,
               showRegionName = input$showRegion,
               type = input$selectTreeType)
    })
  })
  observeEvent(input$msa,{
    # add selected regions
    selectedSet <- getSequencesByRegions(input$selectRegion)
    stringSet$data <- unionDNASets(stringSet$data, selectedSet)
    if (length(stringSet$data) < 3) {
      updateActionButton(session, "msa", label = "Comparison fails. Please
                         Provide at least 3 genomes")
    } else {
      shinybusy::show_spinner() # show the spinner
      alignment$data <- multipleSeqAlign(stringSet$data, input$selectAlgorithm)
      shinybusy::hide_spinner() # hide the spinner
      alignment$tree <- createTree(alignment$data)
      plotTwoGraphs()
    }
  })
}

shiny::shinyApp(ui = ui, server = server)
