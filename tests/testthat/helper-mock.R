# tests/testthat/helper-mock.R
# Builds a realistic TMESplitResult for testing without basilisk.

make_mock_result <- function(
    n_celltypes = 15L, k_shared = 2L, k_spec_A = 1L, k_spec_B = 0L,
    n_A = 20L, n_B = 20L,
    seed = NULL
) {
    if (!is.null(seed)) set.seed(seed)

    ct_names <- paste0("ct", sprintf("%02d", seq_len(n_celltypes)))
    groups <- c("A", "B")

    W_shared <- matrix(runif(n_celltypes * k_shared), nrow = n_celltypes,
                       ncol = k_shared,
                       dimnames = list(ct_names,
                                       if (k_shared > 0L) paste0("shared_", seq_len(k_shared)) else character(0)))

    W_spec_A <- if (k_spec_A > 0L) {
        matrix(runif(n_celltypes * k_spec_A), nrow = n_celltypes,
               ncol = k_spec_A,
               dimnames = list(ct_names, paste0("A_", seq_len(k_spec_A))))
    } else {
        matrix(numeric(0), nrow = n_celltypes, ncol = 0L,
               dimnames = list(ct_names, character(0)))
    }
    W_spec_B <- if (k_spec_B > 0L) {
        matrix(runif(n_celltypes * k_spec_B), nrow = n_celltypes,
               ncol = k_spec_B,
               dimnames = list(ct_names, paste0("B_", seq_len(k_spec_B))))
    } else {
        matrix(numeric(0), nrow = n_celltypes, ncol = 0L,
               dimnames = list(ct_names, character(0)))
    }

    k_total <- k_shared + k_spec_A + k_spec_B
    n_total <- n_A + n_B
    H_raw <- matrix(runif(n_total * k_total), nrow = n_total, ncol = k_total)
    H_frac <- H_raw / rowSums(H_raw)
    prog_names <- c(
        if (k_shared > 0L) paste0("shared_", seq_len(k_shared)) else character(0),
        if (k_spec_A > 0) paste0("A_", seq_len(k_spec_A)),
        if (k_spec_B > 0) paste0("B_", seq_len(k_spec_B))
    )
    pat_names <- c(if (n_A > 0L) paste0("A_p", seq_len(n_A)) else character(0), if (n_B > 0L) paste0("B_p", seq_len(n_B)) else character(0))
    dimnames(H_frac) <- list(pat_names, prog_names)

    group_vec <- rep(groups, c(n_A, n_B))
    names(group_vec) <- pat_names

    prog_rows <- list()
    idx <- 1L
    if (k_spec_A > 0L) {
        for (i in seq_len(k_spec_A)) {
            prog_rows[[idx]] <- data.frame(
                group = "A", program_idx = i - 1L,
                recon_asymmetry = 0.1, recon_asymmetry_p = 0.05,
                ve_ratio = 5.0, ve_ratio_p = 0.01,
                projection_gap = 0.8, projection_gap_p = 0.03,
                combined_p = 0.02, fdr = 0.02,
                stringsAsFactors = FALSE
            )
            idx <- idx + 1L
        }
    }
    if (k_spec_B > 0L) {
        for (i in seq_len(k_spec_B)) {
            prog_rows[[idx]] <- data.frame(
                group = "B", program_idx = i - 1L,
                recon_asymmetry = 0.05, recon_asymmetry_p = 0.2,
                ve_ratio = 2.0, ve_ratio_p = 0.15,
                projection_gap = 0.3, projection_gap_p = 0.25,
                combined_p = 0.18, fdr = 0.18,
                stringsAsFactors = FALSE
            )
            idx <- idx + 1L
        }
    }
    programs_df <- if (length(prog_rows) > 0L) {
        S4Vectors::DataFrame(do.call(rbind, prog_rows))
    } else {
        S4Vectors::DataFrame(
            group = character(0), program_idx = integer(0),
            recon_asymmetry = numeric(0), recon_asymmetry_p = numeric(0),
            ve_ratio = numeric(0), ve_ratio_p = numeric(0),
            projection_gap = numeric(0), projection_gap_p = numeric(0),
            combined_p = numeric(0), fdr = numeric(0)
        )
    }

    new("TMESplitResult",
        W_shared = W_shared,
        W_specific = list(A = W_spec_A, B = W_spec_B),
        H_fractions = H_frac,
        programs = programs_df,
        p_value = 0.03,
        p_value_per_group = c(A = 0.03, B = 1.0),
        groups = groups,
        k_shared = k_shared,
        k_specific = c(A = k_spec_A, B = k_spec_B),
        call = quote(tmesplit(X, group = g)),
        metadata = list(
            n_perm = 100L,
            primary_combiner = "cauchy",
            group_labels = group_vec
        ))
}
