context('non-standard EM')

u <- v <- t(P1pc[,c(1,3,4,6,9,11,12)])

.runThisTest <- Sys.getenv("RunAllTests") == "yes"

if (.runThisTest) {
  test_that('Return-raw = TRUE works', {
    expect_is(LDS_reconstruction(P1annual, u, v, start.year = 1600, num.restarts = 15, return.raw = TRUE), "list")
  })
}
