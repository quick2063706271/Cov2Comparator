% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readData.R
\name{readGenome}
\alias{readGenome}
\title{Read genome data from user}
\usage{
readGenome(fastaFile, nameToRegionsFile = NULL)
}
\arguments{
\item{fastaFile}{A string indicating the path of .fasta file}

\item{nameToRegionsFile}{A string indicating the path of name to region .csv file}
}
\value{
Returns a DNAStringSet with sequences from fastaFile are add.
}
\description{
A function that read genome sequence .fasta file from user by inputting
directories
}
\examples{
# Example 1
# Load MN985325.1.fasta and MT066156.1.fasta file in inst/extdata
fastaPath1 <- system.file("extdata", "MN985325.1.fasta",
  package="Cov2Comparator")
fastaPath2 <- system.file("extdata", "MT066156.1.fasta",
  package="Cov2Comparator")
nameToRegionsFile <- system.file("extdata", "nameToCountry.txt",
  package="Cov2Comparator")
fasta1 <- readGenome(fastaPath1, nameToRegionsFile)

}
\references{
H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2020). Biostrings:
Efficient manipulation of biological strings.
R package version 2.58.0. https://bioconductor.org/packages/Biostrings
}
