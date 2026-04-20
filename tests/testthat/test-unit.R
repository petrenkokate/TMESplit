# tests/testthat/test-unit.R
# Mocked tests — no Python/basilisk needed. Run in CI.

test_that("TMESplitResult class is registered", {
    res <- new("TMESplitResult")
    expect_s4_class(res, "TMESplitResult")
    expect_equal(res@k_shared, 0L)
    expect_identical(res@groups, character(0))
})

test_that("toy_composition fixture loads correctly", {
    f <- system.file("extdata", "toy_composition.rds", package = "TMESplit")
    expect_true(nzchar(f))
    x <- readRDS(f)
    expect_named(x, c("A", "B"))
    expect_equal(dim(x$A), c(20L, 15L))
    expect_equal(dim(x$B), c(20L, 15L))
    expect_true(all(x$A >= 0))
    expect_true(all(x$B >= 0))
})

test_that("accessors return correct types on mock result", {
    res <- make_mock_result()

    expect_s4_class(programs(res), "DataFrame")
    expect_equal(nrow(programs(res)), 1L)

    expect_type(pvalue(res, "global"), "double")
    expect_length(pvalue(res, "global"), 1L)

    pg <- pvalue(res, "group")
    expect_named(pg, c("A", "B"))

    expect_true(is.matrix(wShared(res)))
    expect_equal(nrow(wShared(res)), 15L)
    expect_equal(ncol(wShared(res)), 2L)

    expect_type(wSpecific(res), "list")
    expect_named(wSpecific(res), c("A", "B"))
    expect_equal(ncol(wSpecific(res)$A), 1L)
    expect_equal(ncol(wSpecific(res)$B), 0L)

    expect_type(hFractions(res), "list")
    expect_named(hFractions(res), c("A", "B"))
    expect_equal(nrow(hFractions(res)$A), 20L)
    expect_equal(nrow(hFractions(res)$B), 20L)
    expect_equal(ncol(hFractions(res)$A), 3L)
    expect_equal(ncol(hFractions(res)$B), 2L)
})

test_that("pvalue() level validation works", {
    res <- make_mock_result()
    expect_error(pvalue(res, level = "nonsense"), "arg")
})

test_that("show() prints without error", {
    res <- make_mock_result()
    expect_output(show(res), "TMESplitResult")
    expect_output(show(res), "k_shared")
    expect_output(show(res), "p \\(global\\)")
})

test_that("show() handles empty result", {
    res <- new("TMESplitResult")
    expect_output(show(res), "TMESplitResult")
})

test_that("matrix method validates group length", {
    X <- matrix(runif(60), nrow = 6, ncol = 10)
    expect_error(
        tmesplit(X, group = c("A", "B")),
        "length\\(group\\) must equal nrow"
    )
})

test_that("data.frame method coerces to matrix", {
    df <- data.frame(a = 1:3, b = 4:6)
    expect_error(
        tmesplit(df, group = c("X", "X")),
        "length\\(group\\) must equal nrow"
    )
})

test_that("SE method validates group_col", {
    skip_if_not_installed("SummarizedExperiment")
    mat <- matrix(runif(20), nrow = 4, ncol = 5)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(props = mat),
        colData = S4Vectors::DataFrame(cond = c("A", "A", "B", "B", "B"))
    )
    expect_error(tmesplit(se, group_col = "nonexistent"), "group_col")
})

test_that("mock result with no specific programs", {
    res <- make_mock_result(k_spec_A = 0L, k_spec_B = 0L)
    expect_equal(res@k_specific, c(A = 0L, B = 0L))
    expect_equal(nrow(programs(res)), 0L)
    expect_equal(ncol(hFractions(res)$A), 2L)
    expect_equal(ncol(hFractions(res)$B), 2L)
})
