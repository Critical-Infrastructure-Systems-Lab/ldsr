context('non-standard EM')

u <- v <- t(NPpc)

.runThisTest <- Sys.getenv("RunAllTests") == "yes"

if (.runThisTest) {
  test_that('Return-raw = TRUE works', {
    expect_is(LDS_reconstruction(NPannual, u, v, start.year = 1200, num.restarts = 2, return.raw = TRUE), "list")
  })
}
