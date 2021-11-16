library(Cov2Comparator)

test_that("test simple multiple sequence alignment", {
  suppressPackageStartupMessages(library(Biostrings))
  set1 <- Biostrings::DNAStringSet("ATCGATCG")
  set2 <- Biostrings::DNAStringSet("ATTTTTTT")
  set3 <- Biostrings::DNAStringSet("ATCGATTT")
  set <- union(set1, set2)
  set <- union(set, set3)
  align <- Biostrings::unmasked(multipleSeqAlign(set))
  expect_length(align, 3)
})

test_that("test getSequenceByRegion using Muscle", {
  suppressPackageStartupMessages(library(Biostrings))
  set1 <- Biostrings::DNAStringSet("ATCGATCG")
  set2 <- Biostrings::DNAStringSet("ATTTTTTT")
  set3 <- Biostrings::DNAStringSet("ATCGATTT")
  set <- union(set1, set2)
  set <- union(set, set3)
  align <- Biostrings::unmasked(multipleSeqAlign(set, algorithm = "Muscle"))
  expect_length(align, 3)
})

test_that("test getSequenceByRegion using wrong sequences", {
  suppressPackageStartupMessages(library(Biostrings))
  set <- "AASSAASS"
  expect_error(multipleSeqAlign(set))
})
