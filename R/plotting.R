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
    H <- x@H_fractions
    group_labels <- x@metadata$group_labels
    if (is.null(group_labels)) {
        groups <- x@groups
        n_per <- nrow(H) %/% length(groups)
        group_labels <- rep(groups, each = n_per)
    }

    df <- data.frame(
        patient = rep(rownames(H), ncol(H)),
        program = rep(colnames(H), each = nrow(H)),
        activity = as.vector(H),
        group = rep(group_labels, ncol(H)),
        stringsAsFactors = FALSE
    )
    df$program <- factor(df$program, levels = colnames(H))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$program,
                                           y = .data$activity,
                                           fill = .data$group)) +
        ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
        ggplot2::labs(x = "Program", y = "Activity (H fraction)",
                      title = "Program activity per patient",
                      fill = "Group") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

    print(p)
    invisible(p)
})

#' @describeIn plotNetwork Method for [TMESplitResult].
#' @param top_k Number of top-loading cell types per program to connect.
#' @param threshold Minimum normalized loading for an edge.
#' @export
setMethod("plotNetwork", "TMESplitResult", function(x, top_k = 3L,
                                                     threshold = 0.3, ...) {
    W_parts <- list()
    if (ncol(x@W_shared) > 0L) W_parts[[1L]] <- x@W_shared
    for (Wg in x@W_specific) {
        if (ncol(Wg) > 0L) W_parts[[length(W_parts) + 1L]] <- Wg
    }
    if (length(W_parts) == 0L)
        stop("No programs to plot.", call. = FALSE)

    W <- do.call(cbind, W_parts)
    Wn <- sweep(W, 2, sqrt(colSums(W^2)) + 1e-12, "/")

    ct_names <- rownames(W)
    n <- length(ct_names)
    top_k <- min(top_k, n)

    edges <- data.frame(a = integer(0), b = integer(0),
                        weight = numeric(0))
    for (j in seq_len(ncol(Wn))) {
        top_idx <- order(Wn[, j], decreasing = TRUE)[seq_len(top_k)]
        for (ii in seq_along(top_idx)) {
            for (jj in seq_along(top_idx)) {
                a <- top_idx[ii]
                b <- top_idx[jj]
                if (a >= b) next
                if (Wn[a, j] < threshold || Wn[b, j] < threshold) next
                w <- min(Wn[a, j], Wn[b, j])
                existing <- which(edges$a == a & edges$b == b)
                if (length(existing) > 0L) {
                    edges$weight[existing] <- max(edges$weight[existing], w)
                } else {
                    edges <- rbind(edges,
                                   data.frame(a = a, b = b, weight = w))
                }
            }
        }
    }

    angles <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
    node_df <- data.frame(
        x = cos(angles), y = sin(angles),
        label = ct_names,
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot()

    if (nrow(edges) > 0L) {
        edge_df <- data.frame(
            x = node_df$x[edges$a], y = node_df$y[edges$a],
            xend = node_df$x[edges$b], yend = node_df$y[edges$b],
            weight = edges$weight
        )
        p <- p +
            ggplot2::geom_segment(
                data = edge_df,
                ggplot2::aes(x = .data$x, y = .data$y,
                             xend = .data$xend, yend = .data$yend,
                             alpha = .data$weight),
                color = "grey50",
                linewidth = 0.5 + 2.5 * edge_df$weight
            )
    }

    p <- p +
        ggplot2::geom_point(data = node_df,
                            ggplot2::aes(x = .data$x, y = .data$y),
                            size = 4, color = "#e74c3c") +
        ggplot2::geom_text(data = node_df,
                           ggplot2::aes(x = .data$x * 1.12,
                                        y = .data$y * 1.12,
                                        label = .data$label),
                           size = 3) +
        ggplot2::coord_equal() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Cell-type co-loading network") +
        ggplot2::guides(alpha = "none")

    print(p)
    invisible(p)
})

#' @describeIn plotSignificance Method for [TMESplitResult].
#' @export
setMethod("plotSignificance", "TMESplitResult", function(x, ...) {
    labels <- "Global"
    values <- x@p_value
    level <- "Level 1"

    for (g in x@groups) {
        labels <- c(labels, paste0("group: ", g))
        values <- c(values, x@p_value_per_group[[g]])
        level <- c(level, "Level 2")
    }

    prog <- as.data.frame(x@programs)
    if (nrow(prog) > 0L && "combined_p" %in% colnames(prog)) {
        for (i in seq_len(nrow(prog))) {
            row <- prog[i, ]
            labels <- c(labels,
                        paste0(row$group, "/p", row$program_idx))
            values <- c(values, row$combined_p)
            level <- c(level, "Level 3")
        }
    }

    nlp <- -log10(pmax(values, 1e-10))

    df <- data.frame(
        label = factor(labels, levels = labels),
        nlp = nlp,
        level = factor(level, levels = c("Level 1", "Level 2", "Level 3")),
        stringsAsFactors = FALSE
    )

    level_colors <- c("Level 1" = "#2c3e50",
                      "Level 2" = "#3498db",
                      "Level 3" = "#e67e22")

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$label,
                                           y = .data$nlp,
                                           fill = .data$level)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_hline(yintercept = -log10(0.05),
                            color = "red", linetype = "dashed",
                            linewidth = 0.8) +
        ggplot2::scale_fill_manual(values = level_colors) +
        ggplot2::labs(x = NULL,
                      y = expression(-log[10](p)),
                      title = "Hierarchical significance",
                      fill = "Level") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

    print(p)
    invisible(p)
})

.plot_not_implemented <- function(what) {
    stop(
        what, "() is scaffolded but not yet implemented. ",
        "Plotting lands in Phase 3.2c; see task_plan.md in the ",
        "EcoSplit paper repo.",
        call. = FALSE
    )
}
