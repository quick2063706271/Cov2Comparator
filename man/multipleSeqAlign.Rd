% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequenceComparison.R
\name{multipleSeqAlign}
\alias{multipleSeqAlign}
\title{Compare multiple sequence alignments}
\usage{
multipleSeqAlign(sequences, algorithm = "clustalw")
}
\arguments{
\item{sequences}{A DNAStringset containing several sequences}

\item{algorithm}{A string indicating the algorithm that user wants to use
to calculate multiple sequence alignment}
}
\value{
Returns a MsaAAMultipleAlignment object which is the result of
multiple sequence alignment
}
\description{
A function that calculates multiple sequence alignments
}
\examples{
# Example 1
# Create a basic msa and then plot a tree from it
library(Biostrings)
set1 <- Biostrings::DNAStringSet("ATCGATCG")
set2 <- Biostrings::DNAStringSet("ATTTTTTT")
set3 <- Biostrings::DNAStringSet("ATCGATTT")
setTotal <- union(set1, set2)
setTotal <- union(setTotal, set3)
align <- multipleSeqAlign(setTotal)
}
\references{
Charif D, Lobry J. 2007. “SeqinR 1.0-2: a contributed package to the R
project for statistical computing devoted to biological sequences retrieval
and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.),
Structural approaches to sequence evolution: Molecules, networks,
populations, series Biological and Medical Physics, Biomedical Engineering,
207-232. Springer Verlag, New York.

Paradis E. & Schliep K. 2019. ape 5.0: an environment for modern
phylogenetics and evolutionaryanalyses in R. Bioinformatics 35: 526-528.

U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter (2015) msa:
an R package for multiple sequence alignment. Bioinformatics 31(24):3997-
9999. DOI: 10.1093/bioinformatics/btv176.
}
