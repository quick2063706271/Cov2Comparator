#' Compare multiple sequence alignments
#'
#' A function that calculates multiple sequence alignments
#'
#' @param sequences A DNAStringset containing several sequences
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
#' set1 <- Biostrings::DNAStringSet("ATCGATCG")
#' set2 <- Biostrings::DNAStringSet("ATTTTTTT")
#' set3 <- Biostrings::DNAStringSet("ATCGATTT")
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
#' @importFrom msa msaClustalW
#' @importFrom msa msaMuscle
#'
multiplSeqAlign <- function(sequences, algorithm = "ClustalW") {
  if (class(sequences) != "DNAStringSet") {
    stop("Please provide a DNAStringSet as sequences")
  }
  avaialbleAlgorithm <- c("ClustalW", "Muscle")
  if (!is.element(algorithm, avaialbleAlgorithm)) {
    stop("Please input a valid algorithm from (ClustalW, ClustalOmega, Muscle)")
  }
  alignment <- NULL
  if (algorithm == "Muscle") {
    alignment <- msa::msaMuscle(sequences)
  } else {
    alignment <- msa::msaClustalW(sequences)
  }
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
#' set1 <- Biostrings::DNAStringSet("ATCGATCG")
#' set2 <- Biostrings::DNAStringSet("ATTTTTTT")
#' set3 <- Biostrings::DNAStringSet("ATCGATTT")
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
  if (class(alignment) != 'MsaDNAMultipleAlignment') {
    stop("Please provide a MsaDNAMultipleAlignment object as input")
  }
  if (class(outputName) != 'character') {
    stop("Please input a valid output name")
  }
  saveLocation <- paste(getwd(), "/", sep = "")
  Biostrings::writeXStringSet(as(unmasked(alignment), "XStringSet"),
                  file = paste(saveLocation, outputName, sep = ""))
  return(paste(saveLocation, outputName, sep = ""))
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
#' @importFrom pheatmap pheatmap
#' @importFrom Biostrings unmasked

plotAlignment <- function(alignment, refid, startIdx, endIdx) {
  if (class(alignment) != 'MsaDNAMultipleAlignment') {
    stop("Please provide a MsaDNAMultipleAlignment object as alignment")
  }
  if (class(refid) != "character") {
    stop("Please provide a character object as refid")
  }
  if (class(startIdx) != "numeric" | class(endIdx) != "numeric") {
    stop("Please provide numeric value for startIdx and endIdx")
  }
  if (round(startIdx) != startIdx | round(endIdx) != endIdx) {
    stop("Please provide intgers for startIdx and endIdx")
  }
  if (startIdx >= endIdx) {
    stop("Please provide startIdx that is smaller than endIdx")
  }
  xalign <- Biostrings::unmasked(alignment)
  lengthAlign <- as.integer(nchar(as.character(xalign[1])))
  if (startIdx > lengthAlign | endIdx > lengthAlign) {
    stop("Index out of bounds. Please choose provide startIdx, endIndx that is smaller than length of alignmet")
  }
  refSequence <- xalign[names(xalign) == refid]
  seqMatrix <- sapply(1:length(xalign),function(i){
    as.numeric(as.matrix(xalign[i]) == as.matrix(refSequence))
  })
  seqMatrix <- t(seqMatrix)
  row.names(seqMatrix) <- names(xalign)
  title <- paste("Binary heat map of MSA respect to ", refid, " from ", as.character(startIdx), " to ", as.character(endIdx), sep = "")
  heatmap <- pheatmap::pheatmap(seqMatrix[nrow(seqMatrix):1,startIdx:endIdx], cluster_rows=FALSE,cluster_cols=FALSE, main = title, labels_col = "nucleotide")
  return(heatmap)
}
# [END]
