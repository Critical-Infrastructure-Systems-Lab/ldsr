context('non-standard EM')

u <- v <- t(NPpc[601:813])

.runThisTest <- Sys.getenv("RunAllTests") == "yes"

if (.runThisTest) {
  test_that('Return-raw = TRUE works', {
    expect_is(LDS_reconstruction(NPannual, u, v, start.year = 1800, num.restarts = 2, return.raw = TRUE), "list")
  })
}
