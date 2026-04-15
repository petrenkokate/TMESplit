#' @importFrom methods setClass setGeneric setMethod new slotNames validObject
#'   is as
#' @importFrom utils packageVersion
#' @importFrom basilisk BasiliskEnvironment basiliskStart basiliskStop
#'   basiliskRun
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#'   colData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
NULL

.onLoad <- function(libname, pkgname) {
    invisible(NULL)
}
