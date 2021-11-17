#' Read genome data from user
#'
#' A function that read genome sequence .fasta file from user by inputting
#' directories
#'
#' @param fastaFile A string indicating the path of .fasta file
#' @param nameToRegionsFile A string indicating the path of name to region .csv file
#'
#' @return Returns a DNAStringSet with sequences from fastaFile are add.
#'
#' @examples
#' # Example 1
#' # Load MN985325.1.fasta and MT066156.1.fasta file in inst/extdata
#' fastaPath1 <- system.file("extdata", "MN985325.1.fasta",
#'   package="Cov2Comparator")
#' fastaPath2 <- system.file("extdata", "MT066156.1.fasta",
#'   package="Cov2Comparator")
#' nameToRegionsFile <- system.file("extdata", "nameToCountry.txt",
#'   package="Cov2Comparator")
#' fasta1 <- readGenome(fastaPath1, nameToRegionsFile)
#'
#'
#'
#' @references
#'H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#'Efficient manipulation of biological strings.
#'R package version 2.58.0. https://bioconductor.org/packages/Biostrings
#'
#' @export
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings DNAStringSet
readGenome <- function(fastaFile, nameToRegionsFile = NULL) {
  # check input types and if file exits
  if (class(fastaFile) != "character") {
    stop("Fasta file path should be a string")
  }
  if (! file.exists(fastaFile)) {
    stop("fastaFile not exist")
  }
  userSequence <- Biostrings::readDNAStringSet(fastaFile)
  if (is.null(nameToRegionsFile)) {
    return(userSequence)
  }
  if (! file.exists(nameToRegionsFile)) {
    stop("nameToRegionsFile not exist")
  }
  # read nameToRegionsFile
  nameToRegions <- readNameToRegions(nameToRegionsFile = nameToRegionsFile)
  if (length(userSequence) == 0) {
    stop("Fasta File contains no sequence")
  }
  newNames <- vector(mode='character', length=length(userSequence))
  for (i in seq_along(userSequence)) {
    # assign name as "accessionID region"
    id <- strsplit(names(userSequence)[i], " ")[[1]][1]
    region <- nameToRegions[nameToRegions['Name'] == id][2]
    newNames[i] <- paste(id, region, sep = " ")
  }
  # assign names (ids) to genomes
  names(userSequence) <- newNames
  return(Biostrings::DNAStringSet(userSequence))
}
#' Read a name to region dataframe. A helper function for readGenome
#'
#' A function that read name to region dataframe from user by inputting
#' directories.
#'
#' @param nameToRegionsFile A string indicating the path of name to region .csv file
#'
#' @return Returns a DNAStringSet with sequences from fastaFile are add.
#'
#' @references
#'H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#'Efficient manipulation of biological strings.
#'R package version 2.58.0. https://bioconductor.org/packages/Biostrings
#'
#' @importFrom utils read.csv
#' @keywords internal

readNameToRegions <- function(nameToRegionsFile) {
  # check input
  if (! is.character(nameToRegionsFile)) {
    stop("nameToRegionsFile path should be a string")
  }
  nameToRegions <- utils::read.csv(nameToRegionsFile, header = FALSE)
  names(nameToRegions) <- c("Name", "Region")
  return(nameToRegions)
}


################################################################################


#' Read genome sequence data from user by inputing region
#'
#' A function that takes region as input. It will look up the accessionID of
#' SARS-COV2 of this region. Then it will retrieve the sequence from NCBI
#' database online and return as a DNAStringSet object.
#'
#' @param region A string indicating the interested region in which SARS-COV2
#' are located that you want to study
#'
#' @return Returns a DNAStringSet with sequences.
#'
#' @examples
#' # Example 1
#' # Load MT066156.1 (accessionID of SARS-COV2 in Italy) by
#' # using "Italy" as input
#' region <- "Italy"
#' italySequence <- getSequenceByRegion(region)
#'
#' @references
#'H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#'Efficient manipulation of biological strings.
#'R package version 2.58.0. https://bioconductor.org/packages/Biostrings
#'
#'Paradis E. & Schliep K. (2019). ape 5.0: an environment for modern
#'phylogenetics and evolutionary
#'analyses in R. Bioinformatics 35: 526-528.
#'
#' @export
#' @importFrom Biostrings DNAStringSet
#' @importFrom ape read.GenBank


getSequenceByRegion <- function(region) {
  # check input
  if (! is.character(region)) {
    stop("Please input a character as argument")
  }
  accessionIDToRegion <- get0("accessionIDToRegion",
                              envir = asNamespace("Cov2Comparator"))
  # accessionIDToRegion <- Cov2Comparator:::accessionIDToRegion
  # check if region is valid
  if (! is.element(region, accessionIDToRegion$Region)) {
    stop("Sorry, currently data do not contain your interested region yet. \n
         Please input regions from following: (Canada, Italy, Wuhan, USA,
         Kenya, Bahrain, Germany, Pakistan, Britain,
        Thailand)")
  }
  accessionId <- accessionIDToRegion[accessionIDToRegion['Region'] == region][1]
  # Use ape package to retrieve sequence online
  dnabin <- ape::read.GenBank(accessionId, as.character = TRUE)
  sequenceString <- paste(dnabin[[1]], collapse = "")
  # modify type
  aastringSet <- Biostrings::DNAStringSet(sequenceString)
  names(aastringSet) <- paste(accessionId, region, sep = " ")
  return(Biostrings::DNAStringSet(aastringSet))
}


################################################################################


#' Read multiple genome sequence data from user by inputting a vector of regions
#'
#' A function that takes regions as input. It will look up the accessionIDs of
#' SARS-COV2 of this region. Then it will retrieve the sequences from NCBI
#' database online and return as a DNAStringSet object.
#'
#' @param regions A vector of strings indicating the interested regions in
#' which SARS-COV2 are located that you want to study
#'
#' @return Returns a DNAStringSet with sequences.
#'
#' @examples
#' # Example 1
#' # Load MT066156.1 (accessionID of SARS-COV2 in Italy) by
#' # using "Italy" as input
#' regions <- c("Italy", "Canada", "USA")
#' Sequence <- getSequencesByRegions(regions)
#'
#' @references
#'H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
#'Efficient manipulation of biological strings.
#'R package version 2.58.0. https://bioconductor.org/packages/Biostrings
#'
#'Paradis E. & Schliep K. 2019. ape 5.0: an environment for modern
#'phylogenetics and evolutionary
#'analyses in R. Bioinformatics 35: 526-528.
#'
#' @export
#' @importFrom Biostrings DNAStringSet
getSequencesByRegions <- function(regions) {
  # check input
  if (length(regions) < 1) {
    stop("Please input a vector containing at least one element")
  }
  if (is.element(FALSE, is.element(regions, accessionIDToRegion$Region))) {
    stop("Your input contain region has not been supported yet")
  }
  # get sequence of first region first
  newNames <- vector(mode='character', length=length(regions))
  sequences <- getSequenceByRegion(regions[1])
  newNames[1] <- names(sequences)[1]
  if (length(regions) == 1) {
    return(Biostrings::DNAStringSet(sequences))
  }
  # get sequence for rest of regions
  regions <- regions[-1]
  for (i in seq_along(regions)) {
    iseq <- getSequenceByRegion(regions[i])
    newNames[i + 1] <- names(iseq)[1]
    sequences = union(sequences, iseq)
  }
  names(sequences) <- newNames
  return(Biostrings::DNAStringSet(sequences))
}

# [END]

