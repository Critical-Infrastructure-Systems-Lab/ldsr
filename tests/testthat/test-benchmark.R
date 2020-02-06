context("Benchmark")

test_that("Benchmark works", {

    expect_is(PCR_reconstruction(P1annual, P1pc, start.year = 1600), "list")
})

test_that("Benchmark works with AIC = -Infinity", {

  expect_warning(fit <- PCR_reconstruction(ldsr:::P5annual, ldsr:::P5pc, start.year = 1600))
  expect_is(fit, "list")
})

test_that("Benchmark ensemble works", {
  expect_is(PCR_ensemble(P1annual, list(P1pc, P1pc), start.year = 1600), "list")
})
