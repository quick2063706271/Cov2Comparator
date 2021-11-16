library(Cov2Comparator)

test_that("test getSequenceByRegion using correct country", {
  region <- "Italy"
  italySequence <- getSequenceByRegion(region)
  expect_type(italySequence, "S4")
  expect_s4_class(italySequence, "DNAStringSet")
  expect_length(italySequence, 1)
  expect_equal(names(italySequence), "MT066156.1 Italy")
})

test_that("test getSequenceByRegion using country not included", {
  region1 <- "Korea"
  expect_error(getSequenceByRegion(region1))
  region2 <- 3
  expect_error(getSequenceByRegion(region2))
})

