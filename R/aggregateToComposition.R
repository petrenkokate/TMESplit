#' @describeIn aggregateToComposition Default method for
#'   [SingleCellExperiment::SingleCellExperiment].
#' @export
setMethod(
    "aggregateToComposition", "SingleCellExperiment",
    function(sce, patient_col, celltype_col, ...) {
        stop(
            "aggregateToComposition() is scaffolded but not yet implemented. ",
            "Algorithmic body lands in Phase 3.2b; see task_plan.md in the ",
            "EcoSplit paper repo.",
            call. = FALSE
        )
    }
)
