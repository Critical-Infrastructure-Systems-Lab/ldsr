context("Ensemble model")

u <- v <- t(P1pc[,c(1,3,4,6,9,11,12)])
u.list <- list(u, u)
v.list <- list(v, v)

test_that("Ensemble cross validation works in serial mode", {
  expect_is(cvLDS_ensemble(P1annual, u.list, v.list, start.year = 1600, CV.reps = 2, num.restarts = 2, parallel = FALSE),
            "list")
})

.runThisTest <- Sys.getenv("RunAllTests") == "yes"

if (.runThisTest) {

  test_that("Ensemble model learning works", {
    fit <- LDS_ensemble(P1annual, u.list, v.list, start.year = 1600, num.restarts = 12)
    expect_is(fit, "list")
  })


  test_that("Ensemble cross validation works in parallel", {
    expect_is(cvLDS_ensemble(P1annual, u.list, v.list, start.year = 1600, CV.reps = 2, num.restarts = 2, parallel = TRUE),
              "list")
  })

}
