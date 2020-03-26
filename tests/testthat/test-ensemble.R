context("Ensemble model")

u <- v <- t(P1pc[,c(1,3,4,6,9,11,12)])
u.list <- list(u, u)
v.list <- list(v, v)

test_that("Ensemble model works", {
  expect_is(LDS_reconstruction(P1annual, u.list, v.list, start.year = 1600, num.restarts = 2),
            "list")
})


test_that("Ensemble cross validation works", {
  expect_is(cvLDS(P1annual, u.list, v.list, start.year = 1600, num.restarts = 2,
                  Z = make_Z(P1annual$Qa, nRuns = 2, contiguous = FALSE)),
            "list")
})
