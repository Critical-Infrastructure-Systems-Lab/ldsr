context("PCR")

test_that("PCR works", {
  expect_is(PCR_reconstruction(P1annual, P1pc[, c(1, 3, 4, 6, 9, 11, 12)], start.year = 1600),
            "list")
})

test_that("PCR cross-validation works", {
  expect_is(cvPCR(P1annual, P1pc[, c(1, 3, 4, 6, 9, 11, 12)], start.year = 1600),
            "list")
})

test_that("PCR ensemble works", {
  expect_is(PCR_reconstruction(P1annual,
                               list(P1pc[, c(1, 3, 4, 6, 9, 11, 12)],
                                    P1pc[, c(1, 3, 4, 6, 9, 11, 12)]),
                               start.year = 1600),
            "list")
})

test_that("PCR ensemble cross-validation works", {
  expect_is(cvPCR(P1annual,
                  list(P1pc[, c(1, 3, 4, 6, 9, 11, 12)],
                       P1pc[, c(1, 3, 4, 6, 9, 11, 12)]),
                  start.year = 1600),
            "list")
})
