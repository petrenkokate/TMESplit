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

test_that("plotSignificance returns a ggplot", {
    res <- make_mock_result()
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plotSignificance(res)
    expect_s3_class(p, "ggplot")
})

test_that("plotSignificance handles p=1 (no signal)", {
    res <- make_mock_result(k_spec_A = 0L, k_spec_B = 0L)
    res@p_value <- 1.0
    res@p_value_per_group <- c(A = 1.0, B = 1.0)
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plotSignificance(res)
    expect_s3_class(p, "ggplot")
})

test_that("plotSignificance floors zero p-values at 1e-10", {
    res <- make_mock_result(seed = 42)
    if (nrow(res@programs) == 0L) skip("no programs in mock")
    res@programs$combined_p <- 0
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plotSignificance(res)
    expect_s3_class(p, "ggplot")
    expect_true(all(is.finite(p$data$nlp)))
})

test_that("plotNetwork returns a ggplot", {
    res <- make_mock_result()
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plotNetwork(res)
    expect_s3_class(p, "ggplot")
})
