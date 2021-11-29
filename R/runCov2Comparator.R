#' Launch Shiny App For Package Cov2Comparator
#'
#' A function that launches the shiny app for this package.
#' This shiny app allows user to input sequence .fasta files,
#' or select genome in the menu. Then a multiple sequence
#' alignment and a phylogenetic tree will be calculated.
#' Plots for multiple sequence alignment and phylogenetic tree
#' will be presented.
#' The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value; a shiny page is opened.
#'
#' @examples
#' \dontrun{
#' cov2Comparator::runCov2Comparator()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#' Steipe B., ABC project (.utility 4.07) A Bioinformatics Course: Applied Bioinformatics
#' \href{http://steipe.biochemistry.utoronto.ca/abc/index.php/Bioinformatics_Main_Page}{Link}.
#'
#' @export
#' @importFrom shiny runApp
runCov2Comparator <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "Cov2Comparator")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}

# [END]
