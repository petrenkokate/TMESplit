#' @describeIn aggregateToComposition Method for
#'   [SingleCellExperiment::SingleCellExperiment].
#' @export
setMethod(
    "aggregateToComposition", "SingleCellExperiment",
    function(sce, patient_col, celltype_col, ...) {
        cd <- SummarizedExperiment::colData(sce)
        if (!patient_col %in% colnames(cd))
            stop("patient_col '", patient_col, "' not found in colData",
                 call. = FALSE)
        if (!celltype_col %in% colnames(cd))
            stop("celltype_col '", celltype_col, "' not found in colData",
                 call. = FALSE)

        patients <- as.character(cd[[patient_col]])
        celltypes <- as.character(cd[[celltype_col]])

        tab <- table(patient = patients, celltype = celltypes)
        props <- tab / rowSums(tab)
        mat <- unclass(t(props))
        storage.mode(mat) <- "double"

        unique_patients <- colnames(mat)
        first_idx <- match(unique_patients, patients)
        patient_meta <- cd[first_idx, , drop = FALSE]
        rownames(patient_meta) <- unique_patients
        keep_cols <- setdiff(colnames(patient_meta),
                             c(patient_col, celltype_col))
        patient_meta <- patient_meta[, keep_cols, drop = FALSE]

        SummarizedExperiment::SummarizedExperiment(
            assays = list(proportions = mat),
            colData = patient_meta
        )
    }
)
