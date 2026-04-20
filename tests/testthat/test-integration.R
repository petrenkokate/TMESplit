# tests/testthat/test-integration.R
# Full basilisk roundtrip — skipped in CI, run locally.

skip_on_ci()
skip_if_not_installed("reticulate")

.basilisk_available <- function() {
    tryCatch({
        proc <- basilisk::basiliskStart(TMESplit:::.tmesplit_env)
        on.exit(basilisk::basiliskStop(proc))
        TRUE
    }, error = function(e) FALSE)
}

skip_if_not(.basilisk_available(), "basilisk env not provisioned")

test_that("tmesplit() matrix method returns valid TMESplitResult", {
    toy <- readRDS(
        system.file("extdata", "toy_composition.rds", package = "TMESplit")
    )
    mat <- rbind(toy$A, toy$B)
    group <- rep(c("A", "B"), c(nrow(toy$A), nrow(toy$B)))

    res <- tmesplit(mat, group = group,
                    n_perm = 50L, n_runs = 3L, seed = 42L, verbose = FALSE)

    expect_s4_class(res, "TMESplitResult")
    expect_true(res@k_shared >= 0L)
    expect_equal(length(res@groups), 2L)
    expect_equal(sum(vapply(hFractions(res), nrow, integer(1))), nrow(mat))
    expect_true(res@p_value >= 0 && res@p_value <= 1)
})

test_that("tmesplit() SE method works", {
    toy <- readRDS(
        system.file("extdata", "toy_composition.rds", package = "TMESplit")
    )
    mat <- rbind(toy$A, toy$B)
    group <- rep(c("A", "B"), c(nrow(toy$A), nrow(toy$B)))

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(proportions = t(mat)),
        colData = S4Vectors::DataFrame(condition = group)
    )

    res <- tmesplit(se, group_col = "condition",
                    n_perm = 50L, n_runs = 3L, seed = 42L, verbose = FALSE)

    expect_s4_class(res, "TMESplitResult")
    expect_equal(length(res@groups), 2L)
})

test_that("tmesplit() data.frame method works", {
    toy <- readRDS(
        system.file("extdata", "toy_composition.rds", package = "TMESplit")
    )
    df <- as.data.frame(rbind(toy$A, toy$B))
    group <- rep(c("A", "B"), c(nrow(toy$A), nrow(toy$B)))

    res <- tmesplit(df, group = group,
                    n_perm = 50L, n_runs = 3L, seed = 42L, verbose = FALSE)

    expect_s4_class(res, "TMESplitResult")
})

test_that("matrix and SE methods give same k_shared", {
    toy <- readRDS(
        system.file("extdata", "toy_composition.rds", package = "TMESplit")
    )
    mat <- rbind(toy$A, toy$B)
    group <- rep(c("A", "B"), c(nrow(toy$A), nrow(toy$B)))

    res_mat <- tmesplit(mat, group = group,
                        n_perm = 50L, n_runs = 3L, seed = 42L,
                        verbose = FALSE)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(proportions = t(mat)),
        colData = S4Vectors::DataFrame(condition = group)
    )
    res_se <- tmesplit(se, group_col = "condition",
                       n_perm = 50L, n_runs = 3L, seed = 42L,
                       verbose = FALSE)

    expect_equal(res_mat@k_shared, res_se@k_shared)
    expect_equal(res_mat@k_specific, res_se@k_specific)
    expect_equal(res_mat@p_value, res_se@p_value, tolerance = 1e-6)
})

test_that("all 4 plots render on real result", {
    toy <- readRDS(
        system.file("extdata", "toy_composition.rds", package = "TMESplit")
    )
    mat <- rbind(toy$A, toy$B)
    group <- rep(c("A", "B"), c(nrow(toy$A), nrow(toy$B)))

    res <- tmesplit(mat, group = group,
                    n_perm = 50L, n_runs = 3L, seed = 42L, verbose = FALSE)

    pdf(nullfile())
    on.exit(dev.off())

    ht <- plotPrograms(res)
    expect_s4_class(ht, "Heatmap")

    p1 <- plotActivities(res)
    expect_s3_class(p1, "ggplot")

    p2 <- plotSignificance(res)
    expect_s3_class(p2, "ggplot")

    p3 <- plotNetwork(res)
    expect_s3_class(p3, "ggplot")
})

test_that("numeric parity with Python reference", {
    toy <- readRDS(
        system.file("extdata", "toy_composition.rds", package = "TMESplit")
    )

    proc <- basilisk::basiliskStart(TMESplit:::.tmesplit_env)
    on.exit(basilisk::basiliskStop(proc))

    ref <- basilisk::basiliskRun(proc, function() {
        tme <- reticulate::import("tmesplit")
        np <- reticulate::import("numpy")
        pkl <- reticulate::import("pickle")
        io <- reticulate::import("builtins")

        fixture_path <- system.file(
            "python_fixtures", "toy_micro_reference.pkl",
            package = "TMESplit"
        )
        if (!nzchar(fixture_path)) return(NULL)

        f <- io$open(fixture_path, "rb")
        on.exit(f$close())
        pkl$load(f)
    })

    skip_if(is.null(ref),
            "Python reference fixture not installed")

    mat <- rbind(toy$A, toy$B)
    group <- rep(c("A", "B"), c(nrow(toy$A), nrow(toy$B)))

    res <- tmesplit(mat, group = group,
                    n_perm = 100L, n_runs = 10L, seed = 42L,
                    verbose = FALSE)

    ref_eco <- ref$ecosplit
    expect_equal(res@k_shared, as.integer(ref_eco$k_shared))
})
