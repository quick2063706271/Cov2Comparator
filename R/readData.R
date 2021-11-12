#' Read genome data from user
#'
#' A function that read genome sequence .fasta file from user
#'
#' @param fastaFile A string indicating the path of .fasta file
#' @param regions A vector indicating the region to which each sequence in fastaFile
#'  belong
#' @param genomeDB A dataframe storing sequence and their corresponding region
#'
#' @return Returns a dataframe with sequences from fastaFile are add.
#'
#' @examples
#' # Example 1
#' # Using GeneCounts dataset available with package
#' dim(GeneCounts)
#'
#' # Calculate information criteria value
#' InfCriteriaResults <- InfCriteriaCalculation(loglikelihood = -5080,
#'                                              nClusters = 2,
#'                                              dimensionality = ncol(GeneCounts),
#'                                              observations = nrow(GeneCounts),
#'                                              probability = c(0.5, 0.5))
#' # To obtain BIC value from results
#' InfCriteriaResults$BICresults
#'
#' \dontrun{
#' # Example 2
#' # Obtain an external sample RNAseq dataset
#' library(MBCluster.Seq)
#' data("Count")
#' dim(Count)
#'
#' # Calculate information criteria value
#' InfCriteriaResults <- InfCriteriaCalculation(loglikelihood = -5080,
#'                                              nClusters = 2,
#'                                              dimensionality = ncol(Count),
#'                                              observations = nrow(Count),
#'                                              probability = c(0.5, 0.5))
#' InfCriteriaResults$BICresults
#'}
#' @references
#'Akaike, H. (1973). Information theory and an extension of the maximum
#'likelihood principle. In \emph{Second International Symposium on Information
#'Theory}, New York, NY, USA, pp. 267–281. Springer Verlag. \href{https://link.springer.com/chapter/10.1007/978-1-4612-1694-0_15}{Link}
#'
#'Biernacki, C., G. Celeux, and G. Govaert (2000). Assessing a mixture model for
#'clustering with the integrated classification likelihood. \emph{IEEE Transactions on Pattern
#'Analysis and Machine Intelligence} 22. \href{https://hal.inria.fr/inria-00073163/document}{Link}
#'
#'Schwarz, G. (1978). Estimating the dimension of a model. \emph{The Annals of Statistics} 6, 461–464.
#'\href{https://projecteuclid.org/euclid.aos/1176344136}{Link}.
#'
#'Yaqing, S. (2012). MBCluster.Seq: Model-Based Clustering for RNA-seq
#'Data. R package version 1.0.
#'\href{https://CRAN.R-project.org/package=MBCluster.Seq}{Link}.
#'
#' @export
#' @import seqinr
readGenome <- function(fastaFile, nameToRegionsFile, genomeDB) {
  nameToRegions <- readNameToRegions(nameToRegionsFile = nameToRegionsFile)
  userSequence <- Biostrings::readAAStringSet(fastaFile)
  if (len(userSequence == 0) | len(nameToRegions) == 0 ) {
    # raise error
  }
  for (i in seq_along(1: len(userSequence))) {
    names(userSequence)[i] <- nameToRegions[nameToRegions['Name'] == names(userSequence)[i]][2]
  }
  genomeDB <- union(genomeDB, userSequence)
  return(genomeDB)
}

readNameToRegions <- function(nameToRegionsFile) {
  nameToRegions <- read.csv(nameToRegionsFile, header = FALSE)
  names(nameToRegions) <- c("Name", "Region")
  return(nameToRegions)
}

# [END]

