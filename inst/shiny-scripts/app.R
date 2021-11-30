library(shiny)
ui <- fluidPage(
  # Upload files (fileInput)
  selectInput(inputId = "select",
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

  # Select country/genome (selectInput)
  fileInput(inputId = "fastaFile", label = "Upload your Fasta file",
            accept = c(
              ".fasta"
            )
  ),
  fileInput(inputId = "nameToCountryFile", label = "Upload nameToCountry file",
            accept = c(
              ".txt",
              ".csv"
            )
  ),
  # Output heatmap + plot
)
server <- function(input, output) {
  output$result <- renderText({
    paste("You chose", input$select)
  })
}
shiny::shinyApp(ui = ui, server = server)
