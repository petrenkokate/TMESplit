#' @describeIn tmesplit Method for a [SummarizedExperiment::SummarizedExperiment].
#'   The assay is treated as a `cell_types x patients` composition matrix.
#' @param group_col Column in `colData(object)` giving the group label.
#' @param assay Assay name or index to use (default first).
#' @export
setMethod(
    "tmesplit", "SummarizedExperiment",
    function(object, group_col, assay = 1L, ...) {
        .not_implemented_yet("SummarizedExperiment method")
    }
)

#' @describeIn tmesplit Method for a [SingleCellExperiment::SingleCellExperiment].
#'   First aggregates to a patient x cell-type composition via
#'   [aggregateToComposition], then dispatches to the SE method.
#' @param patient_col Column in `colData(object)` identifying patients.
#' @param celltype_col Column in `colData(object)` identifying cell types.
#' @export
setMethod(
    "tmesplit", "SingleCellExperiment",
    function(object, group_col, patient_col, celltype_col, ...) {
        .not_implemented_yet("SingleCellExperiment method")
    }
)

#' @describeIn tmesplit Method for a raw `matrix`. `group` must be a vector
#'   of length `ncol(object)` giving group labels.
#' @param group Character/factor vector of group labels (matrix/data.frame
#'   dispatch only).
#' @export
setMethod(
    "tmesplit", "matrix",
    function(object, group, ...) {
        .not_implemented_yet("matrix method")
    }
)

#' @describeIn tmesplit Method for a `data.frame`, coerced to a matrix then
#'   dispatched to the matrix method.
#' @export
setMethod(
    "tmesplit", "data.frame",
    function(object, group, ...) {
        .not_implemented_yet("data.frame method")
    }
)

.not_implemented_yet <- function(what) {
    stop(
        "TMESplit ", what, " is scaffolded but not yet implemented. ",
        "This is the Phase 3.2a pre-release skeleton. ",
        "Algorithmic body lands in Phase 3.2b; see task_plan.md in the ",
        "EcoSplit paper repo.",
        call. = FALSE
    )
}
