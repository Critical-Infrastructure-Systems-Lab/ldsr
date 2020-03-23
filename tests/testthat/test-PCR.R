context("Benchmark")

test_that("Benchmark works", {

    expect_is(PCR_reconstruction(P1annual, P1pc[, c(1, 3, 4, 6, 9, 11, 12)], start.year = 1600), "list")
})

test_that("Benchmark cross-validation works", {

  expect_is(cvPCR(P1annual, P1pc[, c(1, 3, 4, 6, 9, 11, 12)], start.year = 1600), "list")
})

test_that("Benchmark ensemble works", {
  expect_is(PCR_ensemble(P1annual, list(P1pc, P1pc), start.year = 1600), "list")
})
