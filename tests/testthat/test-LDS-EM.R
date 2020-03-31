context("EM")

u <- v <- t(P1pc[,c(1,3,4,6,9,11,12)])

test_that("Fixed restarts works", {
  fit <- LDS_reconstruction(P1annual, u, v, start.year = 1600, init = make_init(nrow(u), nrow(v), 1))
  expect_is(fit, "list")
})

test_that("Randomized restarts works", {
  fit <- LDS_reconstruction(P1annual, u, v, start.year = 1600, num.restarts = 2)
  expect_is(fit, "list")
})

test_that("Cross validation works with randomized restarts", {
  cv <- cvLDS(P1annual, u, v, start.year = 1600, num.restarts = 2, Z = make_Z(P1annual$Qa, 2))
  expect_is(cv, "list")
})

test_that("u can be NULL", {
  fit <- LDS_reconstruction(P1annual, u = NULL, v, start.year = 1600, num.restarts = 2)
  expect_is(fit, "list")
  cv <- cvLDS(P1annual, u = NULL, v, start.year = 1600, num.restarts = 2, Z = make_Z(P1annual$Qa, 2))
  expect_is(cv, "list")
})

test_that("v can be NULL", {
  fit <- LDS_reconstruction(P1annual, u, v = NULL, start.year = 1600, num.restarts = 2)
  expect_is(fit, "list")
  cv <- cvLDS(P1annual, u, v = NULL, start.year = 1600, num.restarts = 2, Z = make_Z(P1annual$Qa, 2))
  expect_is(cv, "list")
})

test_that("u and v can have different nrows", {
  fit <- LDS_reconstruction(P1annual, u[1:2, ], v, start.year = 1600, num.restarts = 2)
  expect_is(fit, "list")
  cv <- cvLDS(P1annual, u[1:2, ], v, start.year = 1600, num.restarts = 2, Z = make_Z(P1annual$Qa, 2))
  expect_is(cv, "list")
})
