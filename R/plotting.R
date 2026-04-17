#' @describeIn plotPrograms Method for [TMESplitResult].
#' @export
setMethod("plotPrograms", "TMESplitResult", function(x, ...) {
    W_shared <- x@W_shared
    W_spec <- x@W_specific
    groups <- x@groups

    columns <- list()
    col_labels <- character(0)
    col_type <- character(0)

    if (ncol(W_shared) > 0L) {
        columns[[length(columns) + 1L]] <- W_shared
        labs <- paste0("shared_", seq_len(ncol(W_shared)))
        col_labels <- c(col_labels, labs)
        col_type <- c(col_type, rep("shared", ncol(W_shared)))
    }
    for (g in groups) {
        Wg <- W_spec[[g]]
        if (ncol(Wg) > 0L) {
            columns[[length(columns) + 1L]] <- Wg
            labs <- paste0(g, "_", seq_len(ncol(Wg)))
            col_labels <- c(col_labels, labs)
            col_type <- c(col_type, rep(g, ncol(Wg)))
        }
    }

    if (length(columns) == 0L)
        stop("No programs to plot.", call. = FALSE)

    W <- do.call(cbind, columns)
    colnames(W) <- col_labels

    type_colors <- c(shared = "#2c3e50")
    palette <- c("#3498db", "#e74c3c", "#2ecc71", "#9b59b6",
                 "#f39c12", "#1abc9c")
    for (i in seq_along(groups)) {
        type_colors[groups[i]] <- palette[((i - 1L) %% length(palette)) + 1L]
    }

    ha <- ComplexHeatmap::HeatmapAnnotation(
        type = col_type,
        col = list(type = type_colors[unique(col_type)]),
        show_legend = TRUE,
        show_annotation_name = FALSE
    )

    ht <- ComplexHeatmap::Heatmap(
        W,
        name = "loading",
        col = grDevices::colorRampPalette(
            c("#000004", "#420a68", "#932667", "#dd513a", "#fca50a",
              "#fcffa4"))(100),
        top_annotation = ha,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 8),
        column_names_gp = grid::gpar(fontsize = 9),
        column_title = "TMESplit programs",
        ...
    )

    ComplexHeatmap::draw(ht)
    invisible(ht)
})

#' @describeIn plotActivities Method for [TMESplitResult].
#' @export
setMethod("plotActivities", "TMESplitResult", function(x, ...) {
    .plot_not_implemented("plotActivities")
})

#' @describeIn plotNetwork Method for [TMESplitResult].
#' @export
setMethod("plotNetwork", "TMESplitResult", function(x, ...) {
    .plot_not_implemented("plotNetwork")
})

#' @describeIn plotSignificance Method for [TMESplitResult].
#' @export
setMethod("plotSignificance", "TMESplitResult", function(x, ...) {
    .plot_not_implemented("plotSignificance")
})

.plot_not_implemented <- function(what) {
    stop(
        what, "() is scaffolded but not yet implemented. ",
        "Plotting lands in Phase 3.2c; see task_plan.md in the ",
        "EcoSplit paper repo.",
        call. = FALSE
    )
}
