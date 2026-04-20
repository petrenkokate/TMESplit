# tests/testthat/test-plots.R

# ---------- plotPrograms ----------

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

test_that("plotPrograms honours scale / top_n / transpose / program_names", {
    skip_if_not_installed("ComplexHeatmap")
    res <- make_mock_result()
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)

    ht1 <- plotPrograms(res, scale = FALSE, cluster_rows = FALSE,
                        cluster_columns = FALSE)
    expect_s4_class(ht1, "Heatmap")

    ht2 <- plotPrograms(res, top_n = 5L)
    expect_s4_class(ht2, "Heatmap")

    ht3 <- plotPrograms(res, transpose = TRUE)
    expect_s4_class(ht3, "Heatmap")

    ht4 <- plotPrograms(res,
                        program_names = list(shared = c("S1", "S2"),
                                             A = "Aprog"))
    expect_s4_class(ht4, "Heatmap")
})

# ---------- plotActivities ----------

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

test_that("plotActivities supports all test_method options", {
    res <- make_mock_result(seed = 11L)
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    for (tm in c("wilcoxon", "t", "kruskal")) {
        p <- plotActivities(res, test_method = tm)
        expect_s3_class(p, "ggplot")
    }
})

test_that("plotActivities annotates shared programs with per-program p-values", {
    res <- make_mock_result(seed = 7L)
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plotActivities(res, pvalue_format = "%.2f")
    label_layers <- vapply(p$layers, function(l) {
        inherits(l$geom, "GeomText")
    }, logical(1))
    expect_true(any(label_layers))
})

# ---------- plotSignificance ----------

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

# ---------- plotNetwork ----------

test_that("plotNetwork returns a ggplot", {
    res <- make_mock_result()
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plotNetwork(res)
    expect_s3_class(p, "ggplot")
})

test_that("plotNetwork errors when all programs are empty", {
    res <- make_mock_result(k_shared = 0L, k_spec_A = 0L, k_spec_B = 0L)
    expect_error(plotNetwork(res), "No programs to plot")
})

test_that("plotNetwork runs at a stricter threshold", {
    res <- make_mock_result(seed = 3L)
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    p <- plotNetwork(res, top_k = 3L, threshold = 0.25)
    expect_s3_class(p, "ggplot")
})

# ---------- plotSampleHeatmap (new) ----------

test_that("plotSampleHeatmap returns a Heatmap with defaults", {
    skip_if_not_installed("ComplexHeatmap")
    res <- make_mock_result(seed = 5L)
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    ht <- plotSampleHeatmap(res)
    expect_s4_class(ht, "Heatmap")
})

test_that("plotSampleHeatmap honours annotations + annotation_colors", {
    skip_if_not_installed("ComplexHeatmap")
    res <- make_mock_result(seed = 5L)
    samples <- rownames(res@metadata$composition)
    ann <- data.frame(
        Response = rep(c("R", "NR"), length.out = length(samples)),
        row.names = samples, stringsAsFactors = FALSE
    )
    ann_colors <- list(Response = c(R = "#2C7FB8", NR = "#F768A1"))
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    ht <- plotSampleHeatmap(res, annotations = ann,
                            annotation_colors = ann_colors)
    expect_s4_class(ht, "Heatmap")
})

test_that("plotSampleHeatmap respects clustering opt-in", {
    skip_if_not_installed("ComplexHeatmap")
    res <- make_mock_result(seed = 5L)
    grDevices::pdf(nullfile())
    on.exit(grDevices::dev.off(), add = TRUE)
    ht <- plotSampleHeatmap(res, cluster_samples = TRUE,
                            cluster_cell_types = TRUE, row_split = FALSE)
    expect_s4_class(ht, "Heatmap")
})

test_that("plotSampleHeatmap errors when composition is missing", {
    skip_if_not_installed("ComplexHeatmap")
    res <- make_mock_result(seed = 5L)
    res@metadata$composition <- NULL
    expect_error(plotSampleHeatmap(res), "Composition matrix not found")
})
