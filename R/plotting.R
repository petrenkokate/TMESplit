.assemble_W <- function(x, program_names = NULL) {
    W_shared <- x@W_shared
    W_spec <- x@W_specific
    groups <- x@groups

    columns <- list()
    col_labels <- character(0)
    col_type <- character(0)

    auto_names <- function(prefix, n) paste0(prefix, "_", seq_len(n))
    user_names <- function(key, n) {
        if (is.null(program_names) || is.null(program_names[[key]])) {
            NULL
        } else {
            nm <- as.character(program_names[[key]])
            if (length(nm) != n) {
                stop(sprintf(
                    "program_names$%s has length %d but there are %d programs",
                    key, length(nm), n), call. = FALSE)
            }
            nm
        }
    }

    if (ncol(W_shared) > 0L) {
        columns[[length(columns) + 1L]] <- W_shared
        labs <- user_names("shared", ncol(W_shared))
        if (is.null(labs)) labs <- auto_names("shared", ncol(W_shared))
        col_labels <- c(col_labels, labs)
        col_type <- c(col_type, rep("shared", ncol(W_shared)))
    }
    for (g in groups) {
        Wg <- W_spec[[g]]
        if (is.null(Wg) || ncol(Wg) == 0L) next
        columns[[length(columns) + 1L]] <- Wg
        labs <- user_names(g, ncol(Wg))
        if (is.null(labs)) labs <- auto_names(g, ncol(Wg))
        col_labels <- c(col_labels, labs)
        col_type <- c(col_type, rep(g, ncol(Wg)))
    }

    if (length(columns) == 0L)
        stop("No programs to plot.", call. = FALSE)

    W <- do.call(cbind, columns)
    colnames(W) <- col_labels
    list(W = W, col_type = col_type, groups = groups)
}

.zscale_rows <- function(M) {
    Z <- t(scale(t(M)))
    Z[is.na(Z)] <- 0
    Z
}

#' @describeIn plotPrograms Method for [TMESplitResult].
#' @param scale Logical; z-score each row (cell type) across programs before
#'   plotting. Default `TRUE`.
#' @param top_n Integer or `NULL`; if set, keep the union of each program's
#'   top-`top_n` loading cell types (matches the CoVarNet/6_TIME reference
#'   pattern). Default `NULL` keeps all rows.
#' @param transpose Logical; swap axes so programs become rows and cell types
#'   become columns. Default `FALSE`.
#' @param program_names Named list overriding program column labels. Keys are
#'   `"shared"` and each group name; values are character vectors matching the
#'   number of programs in that block. Unset keys fall back to
#'   `shared_<i>` / `<group>_<i>`.
#' @param cluster_rows,cluster_columns Logical; passed through to the
#'   underlying [ComplexHeatmap::Heatmap]. Both `TRUE` by default.
#' @export
setMethod("plotPrograms", "TMESplitResult",
          function(x,
                   scale = TRUE,
                   top_n = NULL,
                   transpose = FALSE,
                   program_names = NULL,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   ...) {
    parts <- .assemble_W(x, program_names = program_names)
    W <- parts$W
    col_type <- parts$col_type
    groups <- parts$groups

    if (isTRUE(scale)) W <- .zscale_rows(W)

    if (!is.null(top_n)) {
        top_n <- as.integer(top_n)
        if (top_n < 1L) stop("top_n must be a positive integer", call. = FALSE)
        keep <- unique(unlist(lapply(seq_len(ncol(W)), function(j) {
            n <- min(top_n, nrow(W))
            order(W[, j], decreasing = TRUE)[seq_len(n)]
        })))
        W <- W[sort(keep), , drop = FALSE]
    }

    type_levels <- c("shared", groups)
    type_levels <- type_levels[type_levels %in% col_type]
    type_palette <- TMESplit::tme_palette[seq_along(type_levels)]
    names(type_palette) <- type_levels

    col_fun <- circlize::colorRamp2(c(-1.5, 0, 1.5), TMESplit::tme_diverging)

    heat_legend <- list(
        title = if (isTRUE(scale)) "Loading\n(z-score)" else "Loading",
        at = c(-1.5, 0, 1.5),
        labels = c("Low", "0", "High"),
        border = "black"
    )

    if (isTRUE(transpose)) {
        M <- t(W)
        row_type <- col_type
        row_ann <- ComplexHeatmap::rowAnnotation(
            type = row_type,
            col = list(type = type_palette),
            show_legend = TRUE,
            show_annotation_name = FALSE
        )
        top_ann <- NULL
        left_ann <- row_ann
    } else {
        M <- W
        top_ann <- ComplexHeatmap::HeatmapAnnotation(
            type = col_type,
            col = list(type = type_palette),
            show_legend = TRUE,
            show_annotation_name = FALSE
        )
        left_ann <- NULL
    }

    ht <- ComplexHeatmap::Heatmap(
        M,
        name = if (isTRUE(scale)) "Loading\n(z-score)" else "Loading",
        col = col_fun,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = grid::gpar(fontsize = 8),
        column_names_gp = grid::gpar(fontsize = 9),
        column_names_rot = 45,
        border = FALSE,
        rect_gp = grid::gpar(col = "white", lwd = 0.5),
        top_annotation = top_ann,
        left_annotation = left_ann,
        heatmap_legend_param = heat_legend,
        column_title = "TMESplit programs",
        ...
    )

    ComplexHeatmap::draw(ht)
    invisible(ht)
})

.build_activity_df <- function(x) {
    H <- x@H_fractions
    groups <- x@groups
    k_shared <- as.integer(x@k_shared)
    k_spec <- x@k_specific
    if (!is.list(H)) {
        stop("H_fractions must be a named list of matrices (one per group); ",
             "got class ", paste(class(H), collapse = "/"),
             ". Did the fixture pre-date the dict layout?",
             call. = FALSE)
    }

    rows <- list()
    shared_labs <- if (k_shared > 0L)
        paste0("shared_", seq_len(k_shared)) else character(0)

    spec_labs_per <- lapply(groups, function(g) {
        ks <- as.integer(k_spec[[g]])
        if (is.na(ks) || ks == 0L) character(0) else paste0(g, "_", seq_len(ks))
    })
    names(spec_labs_per) <- groups

    level_order <- c(shared_labs, unlist(spec_labs_per, use.names = FALSE))

    for (g in groups) {
        Hg <- H[[g]]
        if (is.null(Hg) || length(Hg) == 0L) next
        ng <- nrow(Hg)
        if (k_shared > 0L) {
            for (j in seq_len(k_shared)) {
                rows[[length(rows) + 1L]] <- data.frame(
                    program = shared_labs[j],
                    group = g,
                    activity = Hg[, j],
                    stringsAsFactors = FALSE
                )
            }
        }
        ks <- as.integer(k_spec[[g]])
        if (!is.na(ks) && ks > 0L) {
            for (j in seq_len(ks)) {
                col <- k_shared + j
                rows[[length(rows) + 1L]] <- data.frame(
                    program = spec_labs_per[[g]][j],
                    group = g,
                    activity = Hg[, col],
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    df <- do.call(rbind, rows)
    df$program <- factor(df$program, levels = level_order)
    df$group <- factor(df$group, levels = groups)
    list(df = df,
         level_order = level_order,
         shared_labs = shared_labs,
         spec_labs_per = spec_labs_per)
}

.lookup_program_pvals <- function(x, spec_labs_per) {
    prog_tbl <- tryCatch(as.data.frame(x@programs),
                          error = function(e) data.frame())
    out <- list()
    if (nrow(prog_tbl) == 0L || !"combined_p" %in% colnames(prog_tbl))
        return(out)
    for (g in names(spec_labs_per)) {
        labs <- spec_labs_per[[g]]
        if (length(labs) == 0L) next
        gi <- which(prog_tbl$group == g)
        gi <- gi[order(gi)]
        n <- min(length(gi), length(labs))
        for (i in seq_len(n)) {
            out[[labs[i]]] <- prog_tbl$combined_p[gi[i]]
        }
    }
    out
}

.shared_program_pvals <- function(df, shared_labs, groups, test_method) {
    out <- list()
    if (length(shared_labs) == 0L) return(out)
    if (length(groups) < 2L) return(out)

    test_fn <- switch(test_method,
        "wilcoxon" = function(values, grp) {
            tryCatch(stats::wilcox.test(values ~ grp, exact = FALSE)$p.value,
                     error = function(e) NA_real_)
        },
        "t" = function(values, grp) {
            tryCatch(stats::t.test(values ~ grp)$p.value,
                     error = function(e) NA_real_)
        },
        "kruskal" = function(values, grp) {
            tryCatch(stats::kruskal.test(values ~ grp)$p.value,
                     error = function(e) NA_real_)
        },
        stop("Unknown test_method: ", test_method, call. = FALSE)
    )

    for (lab in shared_labs) {
        sub <- df[df$program == lab, , drop = FALSE]
        if (nrow(sub) < 2L) next
        grp <- as.character(sub$group)
        n_groups <- length(unique(grp))
        if (n_groups < 2L) next
        if (n_groups > 2L && test_method %in% c("wilcoxon", "t")) {
            fn <- function(values, g) {
                tryCatch(stats::kruskal.test(values ~ g)$p.value,
                         error = function(e) NA_real_)
            }
            out[[lab]] <- fn(sub$activity, grp)
        } else {
            out[[lab]] <- test_fn(sub$activity, as.factor(grp))
        }
    }
    out
}

#' @describeIn plotActivities Method for [TMESplitResult].
#' @param test_method Test used for between-group comparison of *shared*
#'   program activities. One of `"wilcoxon"` (default, G=2),
#'   `"t"` (Welch t-test, G=2), or `"kruskal"` (forced for G>=3 regardless).
#'   Group-specific programs cannot be tested between groups (only one
#'   group carries loadings), so they fall back to the Level 3
#'   `combined_p` from the permutation test and get a `(L3)` marker.
#' @param pvalue_format `sprintf` format string for p-values. Default
#'   `"%.3f"`.
#' @param bracket_fraction Height of the bracket tick as a fraction of
#'   the y-axis range. Default `0.015`.
#' @export
setMethod("plotActivities", "TMESplitResult",
          function(x,
                   test_method = c("wilcoxon", "t", "kruskal"),
                   pvalue_format = "%.3f",
                   bracket_fraction = 0.015,
                   ...) {
    test_method <- match.arg(test_method)
    parts <- .build_activity_df(x)
    df <- parts$df
    level_order <- parts$level_order
    shared_labs <- parts$shared_labs
    spec_labs_per <- parts$spec_labs_per

    groups <- levels(df$group)
    fill_palette <- TMESplit::tme_palette[seq_along(groups)]
    names(fill_palette) <- groups

    shared_pvals <- .shared_program_pvals(df, shared_labs, groups,
                                           test_method)
    spec_pvals <- .lookup_program_pvals(x, spec_labs_per)

    y_max <- max(df$activity, na.rm = TRUE)
    y_min <- min(df$activity, na.rm = TRUE)
    y_range <- y_max - y_min
    tick <- y_range * bracket_fraction
    bracket_y <- y_max + tick * 1.5
    label_y <- bracket_y + tick * 1.5
    headroom <- bracket_y + tick * 4

    bracket_rows <- list()
    label_rows <- list()

    for (lab in shared_labs) {
        p <- shared_pvals[[lab]]
        if (is.null(p) || is.na(p)) next
        prog_idx <- which(level_order == lab)
        bracket_rows[[length(bracket_rows) + 1L]] <- data.frame(
            x = prog_idx - 0.2, xend = prog_idx + 0.2,
            y = bracket_y, yend = bracket_y,
            stringsAsFactors = FALSE
        )
        bracket_rows[[length(bracket_rows) + 1L]] <- data.frame(
            x = prog_idx - 0.2, xend = prog_idx - 0.2,
            y = bracket_y, yend = bracket_y - tick,
            stringsAsFactors = FALSE
        )
        bracket_rows[[length(bracket_rows) + 1L]] <- data.frame(
            x = prog_idx + 0.2, xend = prog_idx + 0.2,
            y = bracket_y, yend = bracket_y - tick,
            stringsAsFactors = FALSE
        )
        label_rows[[length(label_rows) + 1L]] <- data.frame(
            program = factor(lab, levels = level_order),
            x = prog_idx,
            y = label_y,
            label = sprintf(paste0("p=", pvalue_format), p),
            stringsAsFactors = FALSE
        )
    }

    for (g in names(spec_labs_per)) {
        for (lab in spec_labs_per[[g]]) {
            p <- spec_pvals[[lab]]
            if (is.null(p) || is.na(p)) next
            prog_idx <- which(level_order == lab)
            label_rows[[length(label_rows) + 1L]] <- data.frame(
                program = factor(lab, levels = level_order),
                x = prog_idx,
                y = label_y,
                label = sprintf(paste0("p=", pvalue_format, " (L3)"), p),
                stringsAsFactors = FALSE
            )
        }
    }

    bracket_df <- if (length(bracket_rows) > 0L)
        do.call(rbind, bracket_rows) else NULL
    label_df <- if (length(label_rows) > 0L)
        do.call(rbind, label_rows) else NULL

    p_plot <- ggplot2::ggplot(df, ggplot2::aes(x = .data$program,
                                                y = .data$activity,
                                                fill = .data$group)) +
        ggplot2::geom_boxplot(outlier.size = 0.8, alpha = 0.85,
                              position = ggplot2::position_dodge(
                                  preserve = "single"))

    if (!is.null(bracket_df)) {
        p_plot <- p_plot +
            ggplot2::geom_segment(data = bracket_df,
                                  ggplot2::aes(x = .data$x, xend = .data$xend,
                                               y = .data$y, yend = .data$yend),
                                  inherit.aes = FALSE,
                                  linewidth = 0.3, color = "black")
    }
    if (!is.null(label_df)) {
        p_plot <- p_plot +
            ggplot2::geom_text(data = label_df,
                               ggplot2::aes(x = .data$x, y = .data$y,
                                            label = .data$label),
                               inherit.aes = FALSE,
                               size = 2.8, fontface = "plain")
    }

    p_plot <- p_plot +
        ggplot2::scale_fill_manual(values = fill_palette) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(
            mult = c(0.02, 0.08)),
            limits = c(NA, headroom)) +
        ggplot2::labs(x = "Program", y = "Activity (H fraction)",
                      title = "Program activity per sample",
                      fill = "Group") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            panel.grid.major.x = ggplot2::element_blank()
        )

    print(p_plot)
    invisible(p_plot)
})

.assemble_W_network <- function(x) {
    groups <- x@groups
    k_shared <- as.integer(x@k_shared)
    k_spec <- x@k_specific

    parts <- list()
    labs <- character(0)

    if (ncol(x@W_shared) > 0L) {
        parts[[length(parts) + 1L]] <- x@W_shared
        labs <- c(labs, paste0("shared_", seq_len(ncol(x@W_shared))))
    }
    for (g in groups) {
        Wg <- x@W_specific[[g]]
        if (is.null(Wg) || ncol(Wg) == 0L) next
        ks <- as.integer(k_spec[[g]])
        parts[[length(parts) + 1L]] <- Wg
        labs <- c(labs, paste0(g, "_", seq_len(ks)))
    }
    if (length(parts) == 0L) stop("No programs to plot.", call. = FALSE)

    W <- do.call(cbind, parts)
    colnames(W) <- labs
    W
}

#' @describeIn plotNetwork Method for [TMESplitResult]. Draws a
#'   Fruchterman-Reingold force-directed graph of cell types.
#'
#'   **Nodes** are cell types; node color = the cell type's dominant
#'   program (the column of the normalised
#'   `[W_shared | W_specific_by_group]` matrix where it has its highest
#'   loading). **Edges** encode *co-loading on the same program*: for each
#'   program column, the top-`top_k` cell types are enumerated; every pair
#'   whose L2-normalised loadings are both at least `threshold` receives an
#'   edge weighted by `min(loading_a, loading_b)`. Same-program edges are
#'   up-weighted 2x in the layout so clusters separate spatially while
#'   cross-program co-loadings remain visible as bridges.
#'
#'   If real-data programs concentrate on one or two cell types (common in
#'   sparse TME compositions), raising `top_k` rarely helps — lower
#'   `threshold` instead (e.g., `0.05`–`0.15`) to expose secondary
#'   contributors. The default `threshold = 0.1` was tuned for real-data
#'   diffuse programs; simulations with sharper per-program concentration
#'   tolerate `0.2`–`0.3`.
#' @param top_k Number of top-loading cell types per program considered for
#'   edges. Default `5`.
#' @param threshold Minimum L2-normalised loading for a cell type to
#'   participate in an edge from that program. Default `0.1`.
#' @param edge_scale Multiplier for edge thickness (default `3`).
#' @param node_scale Multiplier for node size (default `1`).
#' @param seed Random seed for the layout.
#' @export
setMethod("plotNetwork", "TMESplitResult",
          function(x,
                   top_k = 5L,
                   threshold = 0.1,
                   edge_scale = 3,
                   node_scale = 1,
                   seed = 0L,
                   ...) {
    W <- .assemble_W_network(x)
    Wn <- sweep(W, 2, sqrt(colSums(W^2)) + 1e-12, "/")

    prog_labels <- colnames(Wn)
    ct_names <- rownames(Wn)
    n <- length(ct_names)
    top_k <- min(top_k, n)

    dom_idx <- max.col(Wn, ties.method = "first")
    dom_prog <- prog_labels[dom_idx]

    edge_df <- data.frame(a = integer(0), b = integer(0),
                          weight = numeric(0), same = logical(0))
    for (j in seq_len(ncol(Wn))) {
        top_ix <- order(Wn[, j], decreasing = TRUE)[seq_len(top_k)]
        for (ii in seq_along(top_ix)) {
            for (jj in seq_along(top_ix)) {
                a <- top_ix[ii]; b <- top_ix[jj]
                if (a >= b) next
                if (Wn[a, j] < threshold || Wn[b, j] < threshold) next
                w <- min(Wn[a, j], Wn[b, j])
                same <- dom_prog[a] == dom_prog[b] &&
                    dom_prog[a] == prog_labels[j]
                existing <- which(edge_df$a == a & edge_df$b == b)
                if (length(existing) > 0L) {
                    if (w > edge_df$weight[existing]) {
                        edge_df$weight[existing] <- w
                    }
                    edge_df$same[existing] <-
                        edge_df$same[existing] || same
                } else {
                    edge_df <- rbind(edge_df,
                                     data.frame(a = a, b = b, weight = w,
                                                same = same))
                }
            }
        }
    }

    g <- igraph::make_empty_graph(n = n, directed = FALSE)
    igraph::V(g)$name <- ct_names
    if (nrow(edge_df) > 0L) {
        edges_flat <- as.vector(t(as.matrix(edge_df[, c("a", "b")])))
        g <- igraph::add_edges(g, edges_flat,
                               attr = list(
                                   weight = ifelse(edge_df$same,
                                                   edge_df$weight * 2,
                                                   edge_df$weight),
                                   co_weight = edge_df$weight,
                                   same = edge_df$same
                               ))
    }

    set.seed(seed)
    if (igraph::gsize(g) > 0L) {
        layout <- igraph::layout_with_fr(g,
                                          weights = igraph::E(g)$weight,
                                          niter = 500L)
    } else {
        angles <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
        layout <- cbind(cos(angles), sin(angles))
    }
    layout <- scale(layout, center = TRUE, scale = FALSE)
    max_r <- max(sqrt(rowSums(layout^2)))
    if (is.finite(max_r) && max_r > 0) layout <- layout / max_r

    node_df <- data.frame(
        x = layout[, 1], y = layout[, 2],
        label = ct_names,
        program = factor(dom_prog, levels = unique(dom_prog)),
        stringsAsFactors = FALSE
    )

    prog_levels <- levels(node_df$program)
    prog_palette <- TMESplit::tme_palette[seq_along(prog_levels)]
    names(prog_palette) <- prog_levels

    p <- ggplot2::ggplot()

    if (nrow(edge_df) > 0L) {
        seg_df <- data.frame(
            x = node_df$x[edge_df$a], y = node_df$y[edge_df$a],
            xend = node_df$x[edge_df$b], yend = node_df$y[edge_df$b],
            weight = edge_df$weight,
            same = edge_df$same
        )
        max_w <- max(seg_df$weight)
        p <- p +
            ggplot2::geom_segment(
                data = seg_df,
                ggplot2::aes(x = .data$x, y = .data$y,
                             xend = .data$xend, yend = .data$yend,
                             alpha = .data$weight,
                             color = .data$same),
                linewidth = 0.3 + edge_scale * seg_df$weight / max_w
            ) +
            ggplot2::scale_color_manual(
                values = c(`TRUE` = "#444444", `FALSE` = "grey70"),
                guide = "none")
    }

    p <- p +
        ggplot2::geom_point(data = node_df,
                            ggplot2::aes(x = .data$x, y = .data$y,
                                         fill = .data$program),
                            shape = 21, color = "white", stroke = 0.8,
                            size = 5 * node_scale) +
        ggplot2::geom_text(data = node_df,
                           ggplot2::aes(x = .data$x,
                                        y = .data$y - 0.045,
                                        label = .data$label),
                           vjust = 1, size = 2.6) +
        ggplot2::scale_fill_manual(values = prog_palette,
                                    name = "Dominant program") +
        ggplot2::coord_equal() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Cell-type co-loading network") +
        ggplot2::guides(alpha = "none") +
        ggplot2::theme(legend.position = "right",
                       legend.title = ggplot2::element_text(size = 9),
                       legend.text = ggplot2::element_text(size = 8))

    print(p)
    invisible(p)
})

.row_order_by_program <- function(x, cell_types) {
    W <- .assemble_W_network(x)
    missing_rows <- base::setdiff(cell_types, rownames(W))
    if (length(missing_rows) > 0L) {
        stop("Cell types present in composition but missing from W: ",
             paste(head(missing_rows, 5), collapse = ", "),
             call. = FALSE)
    }
    W <- W[cell_types, , drop = FALSE]
    Wn <- sweep(W, 2, sqrt(colSums(W^2)) + 1e-12, "/")
    dom_idx <- max.col(Wn, ties.method = "first")
    dom_prog <- colnames(Wn)[dom_idx]
    max_loading <- vapply(seq_len(nrow(Wn)), function(i) Wn[i, dom_idx[i]],
                          numeric(1))
    prog_levels <- colnames(Wn)
    ord <- order(factor(dom_prog, levels = prog_levels), -max_loading)
    data.frame(
        cell_type = cell_types[ord],
        program = factor(dom_prog[ord], levels = prog_levels),
        stringsAsFactors = FALSE
    )
}

.assign_dominant_program <- function(x) {
    H <- x@H_fractions
    groups <- x@groups
    k_shared <- as.integer(x@k_shared)
    k_spec <- x@k_specific

    shared_labs <- if (k_shared > 0L)
        paste0("shared_", seq_len(k_shared)) else character(0)

    sample_ids <- character(0)
    prog_ids <- character(0)
    group_ids <- character(0)

    for (g in groups) {
        Hg <- H[[g]]
        if (is.null(Hg) || nrow(Hg) == 0L) next
        ks <- as.integer(k_spec[[g]])
        spec_labs <- if (!is.na(ks) && ks > 0L)
            paste0(g, "_", seq_len(ks)) else character(0)
        prog_labs <- c(shared_labs, spec_labs)
        dom <- prog_labs[max.col(Hg, ties.method = "first")]
        sample_ids <- c(sample_ids, rownames(Hg))
        prog_ids <- c(prog_ids, dom)
        group_ids <- c(group_ids, rep(g, nrow(Hg)))
    }

    prog_order <- c(shared_labs)
    for (g in groups) {
        ks <- as.integer(k_spec[[g]])
        if (!is.na(ks) && ks > 0L)
            prog_order <- c(prog_order, paste0(g, "_", seq_len(ks)))
    }

    data.frame(sample = sample_ids,
               program = factor(prog_ids, levels = prog_order),
               group = factor(group_ids, levels = groups),
               stringsAsFactors = FALSE)
}

#' @describeIn plotSampleHeatmap Method for [TMESplitResult].
#' @param composition Optional `samples x cell_types` numeric matrix. If
#'   `NULL` (default), uses the composition stored in the result metadata.
#' @param annotations Optional `data.frame` of per-sample annotations;
#'   rownames must match sample labels in `H_fractions`. Columns become
#'   top-annotation tracks.
#' @param annotation_colors Optional named list of color vectors for each
#'   annotation column. Keys match columns of `annotations`; each value is a
#'   named vector mapping category -> color.
#' @param scale Logical; z-score each cell type across samples. Default
#'   `TRUE`.
#' @param cluster_samples,cluster_cell_types Logical; Ward-cluster within
#'   each program slice. Defaults both `FALSE` — rows are sorted by
#'   dominant program (argmax of normalized `[W_shared | W_specific_by_group]`)
#'   and within each program by loading descending, producing a diagonal
#'   layout. Set either to `TRUE` to opt into Ward clustering.
#' @param row_split Logical; split rows by each cell type's dominant program
#'   (mirrors the column split on samples). Default `TRUE`. Ignored when
#'   `cluster_cell_types = TRUE`.
#' @param show_sample_names Logical; label each sample on the x-axis.
#'   Default `FALSE` (sample count is usually high).
#' @export
setMethod("plotSampleHeatmap", "TMESplitResult",
          function(x,
                   composition = NULL,
                   annotations = NULL,
                   annotation_colors = NULL,
                   scale = TRUE,
                   cluster_samples = FALSE,
                   cluster_cell_types = FALSE,
                   row_split = TRUE,
                   show_sample_names = FALSE,
                   ...) {
    if (is.null(composition)) {
        composition <- x@metadata$composition
    }
    if (is.null(composition))
        stop("Composition matrix not found. Pass `composition` or refit with ",
             "a version of tmesplit() that stores it in metadata.",
             call. = FALSE)

    dom <- .assign_dominant_program(x)
    common <- base::intersect(rownames(composition), dom$sample)
    if (length(common) == 0L)
        stop("No overlap between composition rownames and sample labels.",
             call. = FALSE)
    composition <- composition[common, , drop = FALSE]
    dom <- dom[match(common, dom$sample), , drop = FALSE]

    row_info <- .row_order_by_program(x, colnames(composition))

    M <- t(composition)
    if (isTRUE(scale)) {
        M <- t(scale(t(M)))
        M[is.na(M)] <- 0
    }

    if (!isTRUE(cluster_cell_types)) {
        M <- M[row_info$cell_type, , drop = FALSE]
    }

    col_fun <- circlize::colorRamp2(c(-1.5, 0, 1.5), TMESplit::tme_diverging)

    groups <- x@groups
    group_palette <- TMESplit::tme_palette[seq_along(groups)]
    names(group_palette) <- groups

    prog_levels <- levels(dom$program)
    prog_palette <- TMESplit::tme_palette[seq_along(prog_levels) + 1L]
    names(prog_palette) <- prog_levels

    ann_df <- data.frame(program = dom$program, group = dom$group,
                         stringsAsFactors = FALSE)
    ann_colors <- list(program = prog_palette, group = group_palette)

    if (!is.null(annotations)) {
        if (is.null(rownames(annotations)))
            stop("`annotations` must have rownames matching sample labels.",
                 call. = FALSE)
        ann_sub <- annotations[match(common, rownames(annotations)), ,
                               drop = FALSE]
        ann_df <- cbind(ann_df, ann_sub)
        if (!is.null(annotation_colors)) {
            for (nm in names(annotation_colors)) {
                if (nm %in% colnames(ann_sub))
                    ann_colors[[nm]] <- annotation_colors[[nm]]
            }
        }
        for (nm in colnames(ann_sub)) {
            if (!nm %in% names(ann_colors)) {
                lv <- unique(as.character(ann_sub[[nm]]))
                lv <- lv[!is.na(lv)]
                pal <- TMESplit::tme_palette[
                    ((seq_along(lv) + 5L - 1L) %% length(TMESplit::tme_palette)) + 1L]
                names(pal) <- lv
                ann_colors[[nm]] <- pal
            }
        }
    }

    top_ann <- ComplexHeatmap::HeatmapAnnotation(
        df = ann_df,
        col = ann_colors,
        annotation_name_gp = grid::gpar(fontsize = 8),
        show_legend = TRUE,
        show_annotation_name = TRUE
    )

    row_split_arg <- NULL
    if (isTRUE(row_split) && !isTRUE(cluster_cell_types)) {
        row_split_arg <- row_info$program
    }

    ht <- ComplexHeatmap::Heatmap(
        M,
        name = if (isTRUE(scale)) "Composition\n(z-score)" else "Composition",
        col = col_fun,
        cluster_rows = cluster_cell_types,
        cluster_columns = cluster_samples,
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        column_split = dom$program,
        row_split = row_split_arg,
        cluster_column_slices = FALSE,
        cluster_row_slices = FALSE,
        show_row_names = TRUE,
        show_column_names = show_sample_names,
        row_names_gp = grid::gpar(fontsize = 8),
        column_names_gp = grid::gpar(fontsize = 6),
        column_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
        row_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
        border = FALSE,
        top_annotation = top_ann,
        column_title = "Samples split by dominant program",
        ...
    )

    ComplexHeatmap::draw(ht, merge_legend = TRUE)
    invisible(ht)
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

    level_colors <- c("Level 1" = TMESplit::tme_palette[13],
                      "Level 2" = TMESplit::tme_palette[1],
                      "Level 3" = TMESplit::tme_palette[5])

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
