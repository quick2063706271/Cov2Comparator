#' Compare multiple sequence alignments
#'
#' A function that calculates multiple sequence alignments
#'
#' @param sequences A AAStringset containing several sequences
#' @param algorithm A string indicating the algorithm that user wants to use
#' to calculate multiple sequence alignment
#'
#' @return Returns a msa
#'
#' @examples
#' # Example 1
#' # Create a basic msa and then plot a tree from it
#' package(msa)
#' package(Biostrings)
#' set1 <- Biostrings::AAStringSet("ATCGATCG")
#' set2 <- Biostrings::AAStringSet("ATTTTTTT")
#' set3 <- Biostrings::AAStringSet("ATCGATTT")
#' set <- union(set1, set2)
#' set <- union(set, set3)
#' align <- multiplSeqAlign(set)

#' \dontrun{
#' # Example 2
#'}
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
#' @importFrom msa msa
#'
multiplSeqAlign <- function(sequences, algorithm = "ClustalW") {
  avaialbleAlgorithm <- c("ClustalW", "ClustalOmega", "Muscle")
  if (!is.element(algorithm, avaialbleAlgorithm)) {
    stop("Please input a valid algorithm from (ClustalW, ClustalOmega, Muscle)")
  }
  alignment <- msa::msa(sequences, algorithm)
  return(alignment)
}
#' Save alignment to a fasta file
#'
#' A helper function that saves alignment to a fasta file
#'
#' @param alignment A MsaAAMultipleAlignment obtained from multiple sequence
#' alignment
#' @param outputName A string indicating the path that it saves to
#'
#' @return Returns NULL
#'
#' @examples
#' # Example 1
#' # Create a basic msa and then plot a tree from it
#' package(msa)
#' package(Biostrings)
#' set1 <- Biostrings::AAStringSet("ATCGATCG")
#' set2 <- Biostrings::AAStringSet("ATTTTTTT")
#' set3 <- Biostrings::AAStringSet("ATCGATTT")
#' set <- union(set1, set2)
#' set <- union(set, set3)
#' align <- multiplSeqAlign(set)
#' saveAlignmentToFasta(align, "align.fasta")

#' \dontrun{
#' # Example 2
#'}
#' @references
#'Charif D, Lobry J. 2007. “SeqinR 1.0-2: a contributed package to the R
#'project for statistical computing devoted to biological sequences retrieval
#'and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.),
#'Structural approaches to sequence evolution: Molecules, networks,
#'populations, series Biological and Medical Physics, Biomedical Engineering,
#'207-232. Springer Verlag, New York.

#'U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter (2015) msa:
#'an R package for multiple sequence alignment. Bioinformatics 31(24):3997-
#'9999. DOI: 10.1093/bioinformatics/btv176.
#'
#' @importFrom Biostrings writeXStringSet

saveAlignmentToFasta <- function(alignment, outputName) {
  if (class(alignment) != 'MsaAAMultipleAlignment') {
    stop("Please provide a MsaAAMultipleAlignment object as input")
  }
  if (class(outputName) != 'character') {
    stop("Please input a valid output name")
  }
  Biostrings::writeXStringSet(as(unmasked(alignment), "XStringSet"),
                  file=outputName)
  return(NULL)
}

#' Compare multiple sequence alignments
#'
#' A function that calculates multiple sequence alignments
#'
#' @param sequences A AAStringset containing several sequences
#' @param algorithm A string indicating the algorithm that user wants to use
#' to calculate multiple sequence alignment
#'
#' @return Returns a msa
#'
#' @examples
#' # Example 1
#' # Create a basic msa and then plot a tree from it
#' package(Biostrings)
#' set1 <- Biostrings::AAStringSet("ATCGATCG")
#' set2 <- Biostrings::AAStringSet("ATTTTTTT")
#' set3 <- Biostrings::AAStringSet("ATCGATTT")
#' set <- union(set1, set2)
#' set <- union(set, set3)
#' align <- multiplSeqAlign(set)
#' plotAlignment(align)
#'
#' \dontrun{
#' # Example 2
#'}
#' @references
#'Charif D, Lobry J. 2007. “SeqinR 1.0-2: a contributed package to the R
#'project for statistical computing devoted to biological sequences retrieval
#'and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.),
#'Structural approaches to sequence evolution: Molecules, networks,
#'populations, series Biological and Medical Physics, Biomedical Engineering,
#'207-232. Springer Verlag, New York.
#'
#'Venket Raghavan (2021). seqvisr: Biological Sequence Visualization Functions
#'in R. R package version 0.2.5. https://github.com/vragh/seqvisr
#'
#'U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter (2015) msa:
#'an R package for multiple sequence alignment. Bioinformatics 31(24):3997-
#'9999. DOI: 10.1093/bioinformatics/btv176.
#'
#' @export
#' @importFrom seqvisr msavisr

plotAlignment <- function(alignment, outputName = 'align.fasta') {
  saveAlignmentToFasta(alignment, outputName)
  if (!file.exists(outputName)) {
    stop("File not found")
  }
  seqvisr::msavisr(mymsa = outputName,
          myref = "cc")
  return()
}
# [END]
