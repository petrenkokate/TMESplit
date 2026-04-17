#' @importFrom methods setClass setGeneric setMethod new slotNames validObject
#'   is as show
#' @importFrom utils packageVersion
#' @importFrom basilisk BasiliskEnvironment basiliskStart basiliskStop
#'   basiliskRun
#' @importFrom reticulate import py_to_r r_to_py
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames
#'   colData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_bar geom_segment geom_point
#'   geom_text geom_hline labs theme_minimal theme element_text element_blank
#'   coord_equal scale_fill_manual scale_color_manual .data theme_void guides
#' @importFrom grid unit
#' @importClassesFrom S4Vectors DataFrame
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
NULL

.onLoad <- function(libname, pkgname) {
    invisible(NULL)
}
