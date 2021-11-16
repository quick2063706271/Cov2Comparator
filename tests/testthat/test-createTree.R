library(Cov2Comparator)

test_that("test createTree using correct input", {
  suppressPackageStartupMessages(library(Biostrings))
  suppressPackageStartupMessages(library(msa))
  set1 <- Biostrings::DNAStringSet("ATCGATCG")
  set2 <- Biostrings::DNAStringSet("ATTTTTTT")
  set3 <- Biostrings::DNAStringSet("ATCGATTT")
  setTotal <- union(set1, set2)
  setTotal <- union(setTotal, set3)
  align <- msa::msaClustalW(Biostrings::DNAStringSet(setTotal))
  t <- createTree(align)
  expect_length(t$tip.label, 3)
  expect_length(t$Nnode, 1)
})

test_that("test createTree using wrong input", {
  align <- 123
  expect_error(createTree(align))
})


