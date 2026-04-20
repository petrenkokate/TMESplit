#' Run TMESplit on cell-type composition data
#'
#' Generic entry point. Methods dispatch on the class of `object`:
#' [SummarizedExperiment::SummarizedExperiment],
#' [SingleCellExperiment::SingleCellExperiment], `matrix`, or `data.frame`.
#'
#' @param object Input composition data.
#' @param ... Further arguments passed to the Python backend (e.g. `beta_loss`,
#'   `n_perm`, `rerank_in_perm`, `combiners`, `random_state`).
#'
#' @return A [TMESplitResult].
#' @examples
#' # Scaffold: bodies raise an informative error until Phase 3.2b.
#' X <- matrix(runif(60), nrow = 6, ncol = 10)
#' tryCatch(
#'     tmesplit(X, group = rep(c("A", "B"), each = 3)),
#'     error = function(e) message(conditionMessage(e))
#' )
#' @export
setGeneric("tmesplit", function(object, ...) standardGeneric("tmesplit"))

#' Aggregate a SingleCellExperiment to a patient x cell-type composition
#'
#' @param sce A [SingleCellExperiment::SingleCellExperiment].
#' @param patient_col Column in `colData(sce)` identifying patients.
#' @param celltype_col Column in `colData(sce)` identifying cell types.
#' @param ... Unused.
#' @return A [SummarizedExperiment::SummarizedExperiment] with composition
#'   proportions in `assay()` and patient metadata in `colData`.
#' @examples
#' # Scaffold: method body lands in Phase 3.2b.
#' showMethods("aggregateToComposition")
#' @export
setGeneric(
    "aggregateToComposition",
    function(sce, patient_col, celltype_col, ...)
        standardGeneric("aggregateToComposition")
)

#' Accessor: per-program statistics
#' @param x A [TMESplitResult].
#' @return A [S4Vectors::DataFrame] with one row per program.
#' @examples
#' res <- new("TMESplitResult")
#' programs(res)
#' @export
setGeneric("programs", function(x) standardGeneric("programs"))

#' Accessor: p-value at a chosen hierarchy level
#' @param x A [TMESplitResult].
#' @param level One of `"global"`, `"group"`, `"program"`.
#' @return Numeric scalar (`"global"`), named numeric vector (`"group"`),
#'   or a [S4Vectors::DataFrame] (`"program"`).
#' @examples
#' res <- new("TMESplitResult")
#' pvalue(res, level = "global")
#' @export
setGeneric("pvalue", function(x, level = "global") standardGeneric("pvalue"))

#' Accessor: shared-program W matrix
#' @param x A [TMESplitResult].
#' @return A numeric matrix of shape `cell_types x k_shared` (possibly empty).
#' @examples
#' res <- new("TMESplitResult")
#' wShared(res)
#' @export
setGeneric("wShared", function(x) standardGeneric("wShared"))

#' Accessor: group-specific W matrices
#' @param x A [TMESplitResult].
#' @return A named `list` of numeric matrices, one per group.
#' @examples
#' res <- new("TMESplitResult")
#' wSpecific(res)
#' @export
setGeneric("wSpecific", function(x) standardGeneric("wSpecific"))

#' Accessor: patient-level program activities (soft assignment)
#' @param x A [TMESplitResult].
#' @return Named list of per-group H matrices
#'   (`n_patients_g x (k_shared + k_specific_g)`); rows sum to 1.
#' @examples
#' res <- new("TMESplitResult")
#' hFractions(res)
#' @export
setGeneric("hFractions", function(x) standardGeneric("hFractions"))

#' Plot program loadings (W) as a heatmap
#' @param x A [TMESplitResult].
#' @param ... Passed to [ComplexHeatmap::Heatmap].
#' @return A [ComplexHeatmap::Heatmap] object (invisibly), drawn as a side
#'   effect.
#' @examples
#' # Scaffold: body lands in Phase 3.2c.
#' res <- new("TMESplitResult")
#' tryCatch(plotPrograms(res), error = function(e) message(conditionMessage(e)))
#' @export
setGeneric("plotPrograms", function(x, ...) standardGeneric("plotPrograms"))

#' Plot patient-level program activities by group
#' @param x A [TMESplitResult].
#' @param ... Extra aesthetics.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' res <- new("TMESplitResult")
#' tryCatch(plotActivities(res), error = function(e) message(conditionMessage(e)))
#' @export
setGeneric("plotActivities",
           function(x, ...) standardGeneric("plotActivities"))

#' Plot program co-occurrence / cell-type network
#' @param x A [TMESplitResult].
#' @param ... Extra aesthetics.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' res <- new("TMESplitResult")
#' tryCatch(plotNetwork(res), error = function(e) message(conditionMessage(e)))
#' @export
setGeneric("plotNetwork",
           function(x, ...) standardGeneric("plotNetwork"))

#' Plot per-sample composition heatmap split by dominant program
#'
#' Samples (rows) x cell types (columns) composition heatmap, z-scored per
#' cell type, with samples grouped by their dominant program (hard
#' assignment from `H_fractions`).
#'
#' @param x A [TMESplitResult].
#' @param ... Extra aesthetics.
#' @return A [ComplexHeatmap::Heatmap] object (invisibly), drawn as a side
#'   effect.
#' @examples
#' res <- new("TMESplitResult")
#' tryCatch(plotSampleHeatmap(res),
#'          error = function(e) message(conditionMessage(e)))
#' @export
setGeneric("plotSampleHeatmap",
           function(x, ...) standardGeneric("plotSampleHeatmap"))

#' Plot significance (permutation null vs observed)
#' @param x A [TMESplitResult].
#' @param ... Extra aesthetics.
#' @return A [ggplot2::ggplot] object.
#' @examples
#' res <- new("TMESplitResult")
#' tryCatch(plotSignificance(res),
#'          error = function(e) message(conditionMessage(e)))
#' @export
setGeneric("plotSignificance",
           function(x, ...) standardGeneric("plotSignificance"))
