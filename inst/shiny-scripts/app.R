library(shiny)
library(shinybusy)
library(shinyalert)

ui <- fluidPage(
  # Use pop up warning message
  useShinyalert(),

  # Defines color of layout
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

      # Selction manual of regions
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
      tags$p(
        "Example file: Cov2Comparator/inst/extdata/MN985325.1.fasta"
      ),
      fileInput(inputId = "nameToRegionFile", label = "Upload nameToRegion file",
                accept = c(
                  ".txt",
                  ".csv"
                )
      ),
      tags$p(
        "Example file: Cov2Comparator/inst/extdata/nameToCountry.txt"
      ),
      actionButton(inputId = "upload", label = "Upload Genome"),
      br(),
      br(),

      # Select algorithm to run MSA
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

      # Options to modify MSA plots
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

      # Options to modify tree plot
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
      shinyalert("Upload Failed! No files uploaded.", type = "error")
    } else {
      uploadedSet <- readGenome(fastaFile = input$fastaFile$datapath,
                                nameToRegionsFile = input$nameToRegionFile$datapath)
      stringSet$data <- uploadedSet
      shinyalert("Upload Success!!!", type = "success")
      updateActionButton(session, "upload", label = "Upload Success!!!")
    }
  })

  alignment <- reactiveValues(data = NULL)

  startIdx <- reactive({
    return(input$startIdx)
  })

  endIdx <- reactive({
    return(input$endIdx)
  })

  output$plotMsa <- renderPlot({
    if (! is.null(alignment$data)) {
      plotAlignment(alignment$data,
                    startIdx = as.numeric(startIdx()),
                    endIdx = as.numeric(endIdx()))
    }
    })

  output$plotTree <- renderPlot({
    if (! is.null(alignment$data)) {
      plotTree(tree = alignment$tree,
               showRegionName = input$showRegion,
               type = input$selectTreeType)
    }
  })

  observeEvent(input$msa,{
    # add selected regions
    selectedSet <- getSequencesByRegions(input$selectRegion)
    stringSet$data <- unionDNASets(stringSet$data, selectedSet)
    if (length(stringSet$data) < 3) {
      shinyalert("Comparison fails. Please
                         Provide at least 3 genomes", type = "error")
    } else {
      shinyalert("Comparison Success!!! Please wait for few minutes.", type = "success")
      shinybusy::show_spinner() # show the spinner
      alignment$data <- multipleSeqAlign(stringSet$data, input$selectAlgorithm)
      shinybusy::hide_spinner() # hide the spinner
      alignment$tree <- createTree(alignment$data)
      updateActionButton(session, "upload", label = "Upload Genome")
    }
  })
}

shiny::shinyApp(ui = ui, server = server)
