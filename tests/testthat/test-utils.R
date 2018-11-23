context("Utilities")

test_that("Time-shifting works", {

    dt <- data.table(year = rep(2001:2003, each = 12),
                     month = rep(1:12, 3))

    dt1 <- water_to_calendar_year(dt, -3)
    ans1 <- data.table(year = rep(2001:2002, each = 12),
                       month = rep(1:12, 2))

    dt2 <- water_to_calendar_year(dt, 3)
    ans2 <- data.table(year = rep(2002:2003, each = 12),
                       month = rep(1:12, 2))

    expect_equal(dt1, ans1)
    expect_equal(dt2, ans2)
})
