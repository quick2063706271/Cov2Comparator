library(Cov2Comparator)

test_that("test getSequencesByRegions using one country", {
  regions <- c("Italy")
  italySequence <- getSequencesByRegions(regions)
  expect_type(italySequence, "S4")
  expect_s4_class(italySequence, "DNAStringSet")
  expect_length(italySequence, 1)
  expect_equal(names(italySequence), "MT066156.1 Italy")
})

test_that("test getSequencesByRegions using multiple country", {
  regions <- c("Canada", "Italy", "Wuhan", "USA", "Kenya", "Bahrain", "Germany", "Pakistan", "Britain",
                "Thailand")
  seqs <- getSequencesByRegions(regions)
  expect_length(seqs, length(regions))
})

test_that("test getSequencesByRegions using multiple country with not supported countries", {
  regions <- c("Canada", "Italy", "Wuhan", "korea")
  expect_error(getSequencesByRegions(regions))
})
