library(Cov2Comparator)

test_that("test plotAlignment using correct input", {

  set1 <- Biostrings::DNAStringSet("ATCGATCG")
  set2 <- Biostrings::DNAStringSet("ATTTTTTT")
  set3 <- Biostrings::DNAStringSet("ATCGATTT")
  setTotal <- union(set1, set2)
  setTotal <- union(setTotal, set3)
  names(setTotal) <- c("a region1", "b region2", "c region3")
  align <- multipleSeqAlign(Biostrings::DNAStringSet(setTotal))
  refid <- "b region2"
  startIdx <- 1
  endIdx <- 8
  plotResult <- plotAlignment(align, refid, startIdx, endIdx)
  expect_type(plot, "closure")
})

test_that("test createTree using wrong input", {
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(msa))
  set1 <- Biostrings::DNAStringSet("ATCGATCG")
  set2 <- Biostrings::DNAStringSet("ATTTTTTT")
  set3 <- Biostrings::DNAStringSet("ATCGATTT")
  setTotal <- union(set1, set2)
  setTotal <- union(setTotal, set3)
  names(setTotal) <- c("a", "b", "c")
  align <- msa::msaClustalW(Biostrings::DNAStringSet(setTotal))
  refid <- "b"
  startIdx <- 1
  endIdx <- 9
  expect_error(plotAlignment(align, refid, startIdx, endIdx))
  refid <- "b"
  startIdx <- 18
  endIdx <- 9
  expect_error(plotAlignment(align, refid, startIdx, endIdx))
  refid <- "e"
  startIdx <- 18
  endIdx <- 9
  expect_error(plotAlignment(align, refid, startIdx, endIdx))
  refid <- "a"
  startIdx <- -1
  endIdx <- 9
  expect_error(plotAlignment(align, refid, startIdx, endIdx))
})
