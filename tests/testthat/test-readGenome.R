library(Cov2Comparator)

test_that("loading one fasta file", {
  fastaPath1 <- system.file("extdata", "MN985325.1.fasta", package="Cov2Comparator")
  nameToRegionsFile <- system.file("extdata", "nameToCountry.txt", package="Cov2Comparator")
  fasta1 <- readGenome(fastaPath1, nameToRegionsFile)

  expect_type(fasta1, "S4")
  expect_s4_class(fasta1, "DNAStringSet")
  expect_length(fasta1, 1)
  expect_equal(as.integer(nchar(as.character(fasta1))), 29882)
  expect_equal(names(fasta1), "MN985325.1 USA")
})

test_that("loading file without nameToRegionsFile", {
  fastaPath1 <- system.file("extdata", "MN985325.1.fasta", package="Cov2Comparator")
  nameToRegionsFile <- system.file("extdata", "nameToCountry.txt", package="Cov2Comparator")
  fasta1 <- readGenome(fastaPath1)
  expect_type(fasta1, "S4")
  expect_s4_class(fasta1, "DNAStringSet")
  expect_length(fasta1, 1)
  expect_false(names(fasta1) == "MN985325.1 USA")
})

test_that("Wrong input for fastaFile", {
  fastaPath1 <- "funfufn"
  nameToRegionsFile <- system.file("extdata", "nameToCountry.txt", package="Cov2Comparator")
  expect_error(readGenome(fastaPath1, nameToRegionsFile))
  fastaPath2 <- 12
  expect_error(readGenome(fastaPath1, nameToRegionsFile))
})

test_that("Wrong input for nameToRegionsFile", {
  fastaPath1 <- system.file("extdata", "MN985325.1.fasta", package="Cov2Comparator")
  nameToRegionsFile <- 12
  expect_error(readGenome(fastaPath1, nameToRegionsFile))
  nameToRegionsFile <- 'dsdsd'
  expect_error(readGenome(fastaPath1, nameToRegionsFile))
})
