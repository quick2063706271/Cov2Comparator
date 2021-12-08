library(Cov2Comparator)

test_that("union 2 DNASets", {
  suppressPackageStartupMessages(library(Biostrings))
  set1 <- Biostrings::DNAStringSet("ATCGATCG")
  set2 <- Biostrings::DNAStringSet("ATTTTTTT")
  setTotal <- unionDNASets(set1, set2)

  expect_length(setTotal, 2)
})

test_that("union 1 DNASet with NULL", {
  suppressPackageStartupMessages(library(Biostrings))
  set1 <- Biostrings::DNAStringSet("ATCGATCG")
  set2 <- NULL
  setTotal <- unionDNASets(set1, set2)

  expect_length(setTotal, 1)
})
