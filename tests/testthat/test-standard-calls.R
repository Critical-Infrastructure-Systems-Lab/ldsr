context("EM")

u <- v <- t(P1pc[,c(1,3,4,6,9,11,12)])

test_that("Fixed restarts works", {
  fit <- LDS_reconstruction(P1annual, u, v, start.year = 1600, method = 'EM',
                            init = make_init(nrow(u), 1))
  expect_is(fit, "list")
})

test_that("Randomized restarts works", {
    expect_is(LDS_reconstruction(P1annual, u, v, start.year = 1600, method = 'EM', num.restarts = 2),
              "list")
})

test_that("Cross validation works with fixed restarts", {
    expect_is(cvLDS(P1annual, u, v, start.year = 1600, method = 'EM', init = make_init(nrow(u), 1)),
              "list")
})

test_that("Cross validation works with randomized restarts", {
    expect_is(cvLDS(P1annual, u, v, start.year = 1600, method = 'EM', num.restarts = 2),
              "list")
})
