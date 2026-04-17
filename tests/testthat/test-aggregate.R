test_that("aggregateToComposition produces correct proportions", {
    skip_if_not_installed("SingleCellExperiment")

    n_cells <- 100L
    patients <- sample(paste0("P", 1:5), n_cells, replace = TRUE)
    celltypes <- sample(c("T", "B", "Mono"), n_cells, replace = TRUE)
    counts <- matrix(rpois(n_cells * 10, 5), nrow = 10)

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            patient = patients,
            celltype = celltypes,
            group = ifelse(patients %in% c("P1", "P2", "P3"), "A", "B")
        )
    )

    comp <- aggregateToComposition(sce,
                                   patient_col = "patient",
                                   celltype_col = "celltype")

    expect_s4_class(comp, "SummarizedExperiment")
    expect_equal(ncol(comp), length(unique(patients)))
    expect_equal(nrow(comp), length(unique(celltypes)))
    props <- SummarizedExperiment::assay(comp)
    expect_true(all(props >= 0))
    expect_equal(as.numeric(colSums(props)), rep(1.0, ncol(comp)),
                 tolerance = 1e-12)
})

test_that("aggregateToComposition preserves patient metadata", {
    skip_if_not_installed("SingleCellExperiment")

    patients <- rep(c("P1", "P2", "P3", "P4"), each = 25)
    celltypes <- rep(c("T", "B", "Mono", "NK", "DC"), times = 20)
    grp <- ifelse(patients %in% c("P1", "P2"), "Ctrl", "Treat")
    counts <- matrix(rpois(100 * 10, 5), nrow = 10)

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(
            patient = patients,
            celltype = celltypes,
            condition = grp
        )
    )

    comp <- aggregateToComposition(sce,
                                   patient_col = "patient",
                                   celltype_col = "celltype")

    cd <- SummarizedExperiment::colData(comp)
    expect_true("condition" %in% colnames(cd))
    expect_setequal(rownames(cd), c("P1", "P2", "P3", "P4"))
})

test_that("aggregateToComposition errors on missing columns", {
    skip_if_not_installed("SingleCellExperiment")

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1, nrow = 1, ncol = 5)),
        colData = S4Vectors::DataFrame(x = 1:5)
    )

    expect_error(
        aggregateToComposition(sce, patient_col = "patient",
                               celltype_col = "celltype"),
        "patient_col"
    )
})
