context("EM")

u <- v <- t(P1pc[,c(1,3,4,6,9,11,12)])

test_that("Fixed restarts works and gives a warning for serial computation", {
    expect_warning(fit <- LDS_reconstruction(P1annual, u, v, method = 'EM',
                                             init = c(0.4, rep(-0.03, 7), 0.2, rep(-0.03, 7))))
    expect_is(fit, "list")
})

test_that("Randomized restarts works serially", {
    expect_is(LDS_reconstruction(P1annual, u, v, method = 'EM',
                                 num.restarts = 2,
                                 parallel = FALSE),
              "list")
})

test_that("Randomized restarts works in parallel with num.restarts > 10", {
    skip_on_cran()
    expect_is(LDS_reconstruction(P1annual, u, v, method = 'EM',
                                 num.restarts = 12,
                                 parallel = TRUE),
              "list")
})

test_that("Cross validation works in serial with fixed restarts", {
    expect_is(cvLDS(P1annual, u, v, method = 'EM',
                    init = c(0.4, rep(-0.03, 7), 0.2, rep(-0.03, 7)),
                    CV.reps = 2, parallel = FALSE),
              "list")
})

test_that("Cross validation works in serial with randomized restarts", {
    expect_is(cvLDS(P1annual, u, v, method = 'EM',
                    num.restarts = 2,
                    CV.reps = 2, parallel = FALSE),
              "list")
})

test_that("Cross validation works in parallel with randomized restarts", {
    skip_on_cran()
    expect_is(cvLDS(P1annual, u, v, method = 'EM',
                    num.restarts = 2,
                    CV.reps = 2, parallel = TRUE),
              "list")
})
