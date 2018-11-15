context("Plots")

test_that("plot_reconstruction works", {
    expect_is(plot_reconstruction(ldsr:::P1model$rec, P1annual, 'inst'), "gg")
})
