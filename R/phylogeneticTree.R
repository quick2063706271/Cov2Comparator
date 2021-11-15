#' Create a neighbor join tree object
#'
#' A function that creates a neighbor join tree from multiple sequence alignment
#'
#' @param alignment A MsaAAMultipleAlignment object for creating phylogenetic
#' tree
#' @return Returns an phylo object which is a phylogenetic tree
#'
#' @examples
#' # Example 1
#' # Using GeneCounts dataset available with package
#' dim(GeneCounts)
#'
#' @references
#'Charif D, Lobry J. 2007. “SeqinR 1.0-2: a contributed package to the R
#'project for statistical computing devoted to biological sequences retrieval
#'and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.),
#'Structural approaches to sequence evolution: Molecules, networks,
#'populations, series Biological and Medical Physics, Biomedical Engineering,
#'207-232. Springer Verlag, New York.
#'
#'Paradis E. & Schliep K. 2019. ape 5.0: an environment for modern
#'phylogenetics and evolutionaryanalyses in R. Bioinformatics 35: 526-528.
#'
#'U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter (2015) msa:
#'an R package for multiple sequence alignment. Bioinformatics 31(24):3997-
#'9999. DOI: 10.1093/bioinformatics/btv176.
#'
#' @export
#' @import seqinr ape


createTree <- function(alignment) {
  if (class(alignment) != 'MsaAAMultipleAlignment') {
    stop("Please provide a MsaAAMultipleAlignment object as input")
  }
  hemoAln2 <- msaConvert(alignment, type="seqinr::alignment")
  distanceMatrix <- seqinr::dist.alignment(hemoAln2, "identity")
  as.matrix(distanceMatrix)
  tree <- ape::nj(distanceMatrix)
  return(tree)
}

plotTree <- function(tree, name, showRegionName = TRUE) {
  if (!showRegionName) {
    for (i in seq_along(1: length(hemoTree$tip.label))) {
      tree$tip.label[i] = strsplit(tree$tip.label[[i]], " ")[[1]][1]
    }
  }
  plot(tree, main = name)
  return()
}
# [END]
