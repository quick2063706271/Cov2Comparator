#' Create a neighbor join tree object
#'
#' A function that creates a neighbor join tree from multiple sequence alignment
#'
#' @param alignment A MsaAAMultipleAlignment object for creating phylogenetic
#' tree
#'
#' @return Returns a phylo object which is a phylogenetic tree
#'
#' @examples
#' # Example 1
#' # Create a basic msa and then create a tree from it
#' library(msa)
#' library(Biostrings)
#' set1 <- Biostrings::DNAStringSet("ATCGATCG")
#' set2 <- Biostrings::DNAStringSet("ATTTTTTT")
#' set3 <- Biostrings::DNAStringSet("ATCGATTT")
#' setTotal <- union(set1, set2)
#' setTotal <- union(setTotal, set3)
#' align <- msa::msa(setTotal)
#' createTree(align)
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
#' @importFrom seqinr dist.alignment
#' @importFrom ape nj
#' @importFrom msa msaConvert

createTree <- function(alignment) {
  # check input
  if (class(alignment) != 'MsaDNAMultipleAlignment') {
    stop("Please provide a MsaDNAMultipleAlignment object as input")
  }
  # convert msa to seqinr::alignment
  aln2 <- msa::msaConvert(alignment, type="seqinr::alignment")
  # make a distance matrix from alignment
  distanceMatrix <- seqinr::dist.alignment(aln2, "identity")
  as.matrix(distanceMatrix)
  tree <- ape::nj(distanceMatrix)
  return(tree)
}


################################################################################


#' Plot a neighbor join tree object
#'
#' A function that creates a neighbor join tree from multiple sequence alignment
#'
#' @param tree A phylo object for plotting
#' @param name A char indicating the name of plot
#' @param showRegionName A boolean indicating if showing the region in the plot
#'
#' @return Plot a phylogenetic tree using this phylo object
#'
#' @examples
#' # Example 1
#' # Create a basic msa and then plot a tree from it
#' library(msa)
#' library(Biostrings)
#' set1 <- Biostrings::DNAStringSet("ATCGATCG")
#' set2 <- Biostrings::DNAStringSet("ATTTTTTT")
#' set3 <- Biostrings::DNAStringSet("ATCGATTT")
#' setTotal <- union(set1, set2)
#' setTotal <- union(setTotal, set3)
#' align <- msa::msa(setTotal)
#' tree <- createTree(align)
#' plotTree(tree = tree,
#'          name = "Simple tree",
#'          showRegionName = TRUE)
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
#' @importFrom graphics plot
#'


plotTree <- function(tree, name, showRegionName = TRUE) {
  # check input
  if (class(tree) != "phylo") {
    stop("Not a phylo tree object")
  }
  if (! is.character(name)) {
    stop("Please input a valid plot name")
  }
  if (! showRegionName) {
    for (i in seq_along(1: length(tree$tip.label))) {
      tree$tip.label[i] = strsplit(tree$tip.label[[i]], " ")[[1]][1]
    }
  }
  plot <- graphics::plot(tree, main = name)
  return(plot)
}
# [END]
