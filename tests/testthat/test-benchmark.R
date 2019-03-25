context("Benchmark")

test_that("Benchmark works", {

    expect_is(PCR_reconstruction(P1annual, P1pc), "list")
})

test_that("Benchmark works with AIC = -Infinity", {

  expect_warning(fit <- PCR_reconstruction(ldsr:::P5annual, ldsr:::P5pc))
  expect_is(fit, "list")
})
