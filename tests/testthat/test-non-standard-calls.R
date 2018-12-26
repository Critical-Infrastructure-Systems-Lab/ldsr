context('non-standard EM')

u <- v <- t(P1pc[,c(1,3,4,6,9,11,12)])

test_that('Return-raw = TRUE works', {
    skip_on_cran()
    expect_is(LDS_reconstruction(P1annual, u, v, num.restarts = 15, return.raw = TRUE), "list")
})
