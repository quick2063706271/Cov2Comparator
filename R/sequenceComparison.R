#' Compare multiple sequence alignments
#'
#' A function that calculates multiple sequence alignments
#'
#' @param sequences A DNAStringset containing several sequences
#' @param algorithm A string indicating the algorithm that user wants to use
#' to calculate multiple sequence alignment
#'
#' @return Returns a MsaAAMultipleAlignment object which is the result of
#' multiple sequence alignment
#'
#' @examples
#' # Example 1
#' # Create a basic msa and then plot a tree from it
#' library(Biostrings)
#' set1 <- Biostrings::DNAStringSet("ATCGATCG")
#' set2 <- Biostrings::DNAStringSet("ATTTTTTT")
#' set3 <- Biostrings::DNAStringSet("ATCGATTT")
#' setTotal <- union(set1, set2)
#' setTotal <- union(setTotal, set3)
#' align <- multipleSeqAlign(setTotal)

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
#' @importFrom msa msaClustalW
#' @importFrom msa msaMuscle
#'
multipleSeqAlign <- function(sequences, algorithm = "clustalw") {
  # check input
  if (class(sequences) != "DNAStringSet") {
    stop("Please provide a DNAStringSet as sequences")
  }
  avaialbleAlgorithm <- c("clustalw", "muscle")
  # turn the algorithm to lower case
  lowerAlgorithm <- tolower(algorithm)
  if (! is.element(lowerAlgorithm, avaialbleAlgorithm)) {
    stop("Please input a valid algorithm from (ClustalW, ClustalOmega, Muscle)")
  }
  alignment <- NULL
  if (lowerAlgorithm == "muscle") {
    alignment <- msa::msaMuscle(sequences)
  } else {
    alignment <- msa::msaClustalW(sequences)
  }
  print("==Comparison done!==")
  return(alignment)
}

################################################################################


#' Plot multiple sequence alignments
#'
#' A function that plot a heatmap of multiple sequence alignments, showing all
#' nucleotide that is different from the reference sequence
#'
#' @param alignment A MsaAAMultipleAlignment object which is a result of
#' multiple sequence alignment
#' @param refid A string indicating the name of reference sequence.
#' This sequence must be in the alignment
#' @param startIdx A integer indicating the start position of alignment
#' that user wants to plot
#' @param endIdx A integer indicating the end position of alignment that
#' user wants to plot
#'
#' @return Returns a plot
#'
#' @examples
#' # Example 1
#' # Create a basic msa and then plot a tree from it
#' library(Biostrings)
#' set1 <- Biostrings::DNAStringSet("ATCGATCG")
#' set2 <- Biostrings::DNAStringSet("ATTTTTTT")
#' set3 <- Biostrings::DNAStringSet("ATCGATTT")
#' setTotal <- union(set1, set2)
#' setTotal <- union(setTotal, set3)
#' names(setTotal) <- c("a", "b", "c")
#' align <- multipleSeqAlign(setTotal)
#' plotAlignment(align, refid = "a", startIdx = 1, endIdx = 8)
#'
#' @references
#'Charif D, Lobry J. 2007. “SeqinR 1.0-2: a contributed package to the R
#'project for statistical computing devoted to biological sequences retrieval
#'and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.),
#'Structural approaches to sequence evolution: Molecules, networks,
#'populations, series Biological and Medical Physics, Biomedical Engineering,
#'207-232. Springer Verlag, New York.
#'
#'Raivo Kolde (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12.
#'https://CRAN.R-project.org/package=pheatmap
#'
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom Biostrings unmasked
#' @importFrom BiocGenerics width

plotAlignment <- function(alignment, refid, startIdx, endIdx) {
  # check input
  if (class(alignment) != 'MsaDNAMultipleAlignment') {
    stop("Please provide a MsaDNAMultipleAlignment object as alignment")
  }
  xalign <- Biostrings::unmasked(alignment)
  lengthOfAlign <- BiocGenerics::width(xalign)[1]
  if (missing(startIdx)) {
    startIdx <- 1
  }
  if (missing(endIdx)) {
    endIdx <- lengthOfAlign
  }
  if ((! is.numeric(startIdx)) | (! is.numeric(endIdx))) {
    stop("Please provide numeric value for startIdx and endIdx")
  }
  if (round(startIdx) != startIdx | round(endIdx) != endIdx) {
    stop("Please provide intgers for startIdx and endIdx")
  }
  if (startIdx >= endIdx) {
    stop("Please provide startIdx that is smaller than endIdx")
  }
  if (startIdx < 1 | endIdx < 1) {
    stop("Index out of bounds. Please choose provide startIdx,
         endIndx that is positive number")
  }
  if (startIdx > lengthOfAlign | endIdx > lengthOfAlign) {
    stop("Index out of bounds. Please choose provide startIdx,
         endIndx that is positive number")
  }

  if (is.null(names(xalign))) {
    stop("Your alignment do not have names")
  }
  lengthAlign <- as.integer(nchar(as.character(xalign[1])))
  if (startIdx > lengthAlign | endIdx > lengthAlign) {
    stop("Index out of bounds. Please choose provide startIdx,
         endIndx that is smaller than length of alignmet")
  }
  if (missing(refid)) {
    refid <- names(xalign)[1]
  }
  if (! is.character(class(refid))) {
    stop("Please provide a character object as refid")
  }
  if (! is.element(refid, names(xalign))) {
    stop("Please provide a valid refid")
  }
  # get sequence for reference
  refSequence <- xalign[names(xalign) == refid]
  # make a binary matrix for all genome
  seqMatrix <- sapply(1:length(xalign),function(i){
    as.numeric(as.matrix(xalign[i]) == as.matrix(refSequence))
  })
  # transpose the seqMatrix
  seqMatrix <- t(seqMatrix)
  # assign names to matrix
  row.names(seqMatrix) <- names(xalign)
  title <- paste("Binary heat map of MSA respect to ",
                 refid,
                 " from ",
                 as.character(startIdx),
                 " to ",
                 as.character(endIdx),
                 sep = "")
  heatmap <- pheatmap::pheatmap(seqMatrix[nrow(seqMatrix) : 1,
                                          startIdx : endIdx],
                                cluster_rows=FALSE,
                                cluster_cols=FALSE,
                                main = title,
                                labels_col = "nucleotide")
  return(heatmap)
}
# [END]
