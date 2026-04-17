# tests/testthat/test-plots.R

test_that("plotPrograms returns a Heatmap and draws", {
    skip_if_not_installed("ComplexHeatmap")
    res <- make_mock_result()
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    ht <- plotPrograms(res)
    expect_s4_class(ht, "Heatmap")
})

test_that("plotPrograms works with no specific programs", {
    skip_if_not_installed("ComplexHeatmap")
    res <- make_mock_result(k_spec_A = 0L, k_spec_B = 0L)
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    ht <- plotPrograms(res)
    expect_s4_class(ht, "Heatmap")
})

test_that("plotActivities returns a ggplot and draws", {
    res <- make_mock_result()
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plotActivities(res)
    expect_s3_class(p, "ggplot")
})

test_that("plotActivities works with null-specific result", {
    res <- make_mock_result(k_spec_A = 0L, k_spec_B = 0L)
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plotActivities(res)
    expect_s3_class(p, "ggplot")
})
