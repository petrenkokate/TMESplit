#' @describeIn tmesplit Method for a [SummarizedExperiment::SummarizedExperiment].
#'   The assay should be a `cell_types x patients` composition matrix.
#' @param group_col Column in `colData(object)` giving the group label.
#' @param assay Assay name or index to use (default first).
#' @export
setMethod(
    "tmesplit", "SummarizedExperiment",
    function(object, group_col, assay = 1L, ...) {
        cd <- SummarizedExperiment::colData(object)
        if (!group_col %in% colnames(cd))
            stop("group_col '", group_col, "' not found in colData",
                 call. = FALSE)
        mat <- t(SummarizedExperiment::assay(object, assay))
        group_vec <- as.character(cd[[group_col]])
        .run_tmesplit(mat, group_vec,
                      cell_type_names = rownames(
                          SummarizedExperiment::assay(object, assay)),
                      call_ = match.call(), ...)
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
        comp <- aggregateToComposition(object,
                                       patient_col = patient_col,
                                       celltype_col = celltype_col)
        tmesplit(comp, group_col = group_col, ...)
    }
)

#' @describeIn tmesplit Method for a raw `matrix`. `group` must be a vector
#'   of length `nrow(object)` giving group labels (rows = patients,
#'   cols = cell types).
#' @param group Character/factor vector of group labels (matrix/data.frame
#'   dispatch only).
#' @export
setMethod(
    "tmesplit", "matrix",
    function(object, group, ...) {
        group <- as.character(group)
        if (length(group) != nrow(object))
            stop("length(group) must equal nrow(object)", call. = FALSE)
        ct_names <- colnames(object)
        if (is.null(ct_names)) ct_names <- paste0("ct_", seq_len(ncol(object)))
        .run_tmesplit(object, group, cell_type_names = ct_names,
                      call_ = match.call(), ...)
    }
)

#' @describeIn tmesplit Method for a `data.frame`, coerced to a matrix then
#'   dispatched to the matrix method.
#' @export
setMethod(
    "tmesplit", "data.frame",
    function(object, group, ...) {
        tmesplit(as.matrix(object), group = group, ...)
    }
)

# ---- internal engine --------------------------------------------------------

.run_tmesplit <- function(mat, group_vec, cell_type_names = NULL,
                          call_ = NULL,
                          beta_loss = 0.5,
                          recon_threshold = 0.04,
                          n_perm = 2000L,
                          n_runs = 10L,
                          n_perm_runs = 3L,
                          rerank_in_perm = TRUE,
                          combiners = "cauchy",
                          seed = 42L,
                          n_jobs = 1L,
                          verbose = TRUE,
                          ...) {
    groups <- unique(group_vec)
    freq_dict <- lapply(stats::setNames(groups, groups), function(g) {
        mat[group_vec == g, , drop = FALSE]
    })

    result <- .basilisk_run(function() {
        tme <- reticulate::import("tmesplit")

        freq_py <- lapply(freq_dict, function(m) {
            reticulate::np_array(m, dtype = "float64")
        })

        eco <- tme$ecosplit(
            freq_dict = freq_py,
            cell_type_names = as.list(cell_type_names),
            beta_loss = beta_loss,
            recon_threshold = recon_threshold,
            n_runs = as.integer(n_runs),
            seed = as.integer(seed),
            verbose = verbose,
            ...
        )

        perm <- tme$test_significance(
            freq_dict = freq_py,
            result = eco,
            n_perm = as.integer(n_perm),
            n_runs = as.integer(n_perm_runs),
            beta_loss = beta_loss,
            recon_threshold = recon_threshold,
            combiners = combiners,
            rerank_in_perm = rerank_in_perm,
            n_jobs = as.integer(n_jobs),
            seed = as.integer(seed),
            verbose = verbose
        )

        list(eco = eco, perm = perm)
    })

    .hydrate_result(result$eco, result$perm, groups, group_vec,
                    cell_type_names, call_, beta_loss,
                    freq_dict = freq_dict)
}

.hydrate_result <- function(eco, perm, groups, group_vec,
                            cell_type_names, call_, beta_loss,
                            freq_dict = NULL) {
    k_shared <- as.integer(eco$k_shared)

    k_specific_list <- eco$k_specific
    k_specific <- vapply(groups, function(g) as.integer(k_specific_list[[g]]),
                         integer(1))

    W_shared_raw <- as.matrix(eco$W_shared)
    if (nrow(W_shared_raw) > 0 && ncol(W_shared_raw) > 0) {
        rownames(W_shared_raw) <- cell_type_names
        colnames(W_shared_raw) <- paste0("shared_", seq_len(ncol(W_shared_raw)))
    }

    W_specific <- lapply(stats::setNames(groups, groups), function(g) {
        w <- as.matrix(eco$W_specific[[g]])
        if (nrow(w) > 0 && ncol(w) > 0) {
            rownames(w) <- cell_type_names
            colnames(w) <- paste0(g, "_", seq_len(ncol(w)))
        } else {
            w <- matrix(numeric(0), nrow = length(cell_type_names), ncol = 0L,
                        dimnames = list(cell_type_names, character(0)))
        }
        w
    })

    shared_labs <- if (k_shared > 0L)
        paste0("shared_", seq_len(k_shared)) else character(0)

    H_parts <- list()
    sample_names_per <- list()
    for (g in groups) {
        Hg <- as.matrix(eco$H_fractions[[g]])
        if (!is.null(freq_dict) && !is.null(rownames(freq_dict[[g]]))) {
            sample_labs <- rownames(freq_dict[[g]])
        } else {
            sample_labs <- paste0(g, "_s", seq_len(nrow(Hg)))
        }
        rownames(Hg) <- sample_labs
        ks <- as.integer(k_specific[[g]])
        spec_labs <- if (ks > 0L) paste0(g, "_", seq_len(ks)) else character(0)
        colnames(Hg) <- c(shared_labs, spec_labs)
        H_parts[[g]] <- Hg
        sample_names_per[[g]] <- sample_labs
    }
    H_fractions <- H_parts

    composition <- NULL
    if (!is.null(freq_dict)) {
        comp_parts <- list()
        for (g in groups) {
            Xg <- as.matrix(freq_dict[[g]])
            rownames(Xg) <- sample_names_per[[g]]
            colnames(Xg) <- cell_type_names
            comp_parts[[g]] <- Xg
        }
        composition <- do.call(rbind, comp_parts)
    }

    p_value <- as.numeric(perm$p_value)

    p_per_group_list <- perm$p_value_per_group
    p_value_per_group <- vapply(groups, function(g) {
        as.numeric(p_per_group_list[[g]])
    }, numeric(1))

    prog_table_raw <- perm$programs_table
    if (length(prog_table_raw) > 0L) {
        rows <- lapply(prog_table_raw, function(row) {
            data.frame(
                group = as.character(row$group),
                program_idx = as.integer(row$program_idx),
                recon_asymmetry = as.numeric(row$recon_asymmetry),
                recon_asymmetry_p = as.numeric(row$recon_asymmetry_p),
                ve_ratio = as.numeric(row$ve_ratio),
                ve_ratio_p = as.numeric(row$ve_ratio_p),
                projection_gap = as.numeric(row$projection_gap),
                projection_gap_p = as.numeric(row$projection_gap_p),
                combined_p = as.numeric(row$combined_p),
                fdr = as.numeric(row$fdr),
                stringsAsFactors = FALSE
            )
        })
        programs <- S4Vectors::DataFrame(do.call(rbind, rows))
    } else {
        programs <- S4Vectors::DataFrame(
            group = character(0), program_idx = integer(0),
            recon_asymmetry = numeric(0), recon_asymmetry_p = numeric(0),
            ve_ratio = numeric(0), ve_ratio_p = numeric(0),
            projection_gap = numeric(0), projection_gap_p = numeric(0),
            combined_p = numeric(0), fdr = numeric(0)
        )
    }

    new("TMESplitResult",
        W_shared = W_shared_raw,
        W_specific = W_specific,
        H_fractions = H_fractions,
        programs = programs,
        p_value = p_value,
        p_value_per_group = p_value_per_group,
        groups = groups,
        k_shared = k_shared,
        k_specific = k_specific,
        call = call_,
        metadata = list(
            n_perm = as.integer(perm$n_perm),
            primary_combiner = as.character(perm$primary_combiner),
            beta_loss = beta_loss,
            group_labels = rep(groups,
                               times = vapply(H_parts, nrow, integer(1))),
            sample_names = unlist(sample_names_per, use.names = FALSE),
            composition = composition
        ))
}
