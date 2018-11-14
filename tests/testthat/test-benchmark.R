context("Benchmark")

test_that("Benchmark works", {

    expect_is(PCR_reconstruction(P1annual, P1pc), "list")
})
