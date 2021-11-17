library(Cov2Comparator)

test_that("integration test", {
  set1 <- Biostrings::DNAStringSet("ATCGATCG")
  set2 <- Biostrings::DNAStringSet("ATTTTTTT")
  set3 <- Biostrings::DNAStringSet("ATCGATTT")
  set4 <- getSequenceByRegion("Italy")
  setTotal <- union(set1, set2)
  setTotal <- union(setTotal, set3)
  names(setTotal) <- c("a region1", "b region2", "c region3")
  align <- multipleSeqAlign(Biostrings::DNAStringSet(setTotal))
  tree <- createTree(align)

  expect_length(set4, 1)
  expect_equal(names(set4), "MT066156.1 Italy")
  expect_length(tree$tip.label, 3)
  expect_length(tree$Nnode, 1)
})
