library(Cov2Comparator)

test_that("test plotTree using correct input", {
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(msa))
  set1 <- Biostrings::DNAStringSet("ATCGATCG")
  set2 <- Biostrings::DNAStringSet("ATTTTTTT")
  set3 <- Biostrings::DNAStringSet("ATCGATTT")
  setTotal <- union(set1, set2)
  setTotal <- union(setTotal, set3)
  names(setTotal) <- c("a", "b", "c")
  align <- msa::msa(Biostrings::DNAStringSet(setTotal))
  refid <- "b"
  startIdx <- 1
  endIdx <- 8
  plotResult <- plotAlignment(align, refid, startIdx, endIdx)
  tree <- createTree(align)
  name  <-  "Simple tree"
  showRegionName <- TRUE
  treePlot <- plotTree(tree,
           name,
           showRegionName)
  expect_type(plot, "closure")
  showRegionName <- FALSE
  treePlot <- plotTree(tree,
                       name,
                       showRegionName)
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
  align <- msa::msa(Biostrings::DNAStringSet(setTotal))
  refid <- "b"
  startIdx <- 1
  endIdx <- 8
  plotResult <- plotAlignment(align, refid, startIdx, endIdx)
  tree <- createTree(align)
  name  <-  6
  showRegionName <- TRUE
  expect_error(plotTree(tree,
                       name,
                       showRegionName))
})

