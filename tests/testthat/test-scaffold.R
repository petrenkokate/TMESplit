test_that("TMESplitResult class is registered", {
    res <- new("TMESplitResult")
    expect_s4_class(res, "TMESplitResult")
    expect_equal(res@k_shared, 0L)
    expect_identical(res@groups, character(0))
})

test_that("toy_composition fixture loads", {
    f <- system.file("extdata", "toy_composition.rds", package = "TMESplit")
    expect_true(nzchar(f))
    x <- readRDS(f)
    expect_named(x, c("A", "B"))
    expect_equal(dim(x$A), c(20L, 15L))
    expect_equal(dim(x$B), c(20L, 15L))
})

test_that("tmesplit() scaffold errors with a clear message", {
    X <- matrix(runif(60), nrow = 6, ncol = 10)
    expect_error(
        tmesplit(X, group = rep(c("A", "B"), each = 3)),
        "not yet implemented"
    )
})
