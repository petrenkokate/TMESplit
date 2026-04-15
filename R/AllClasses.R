#' TMESplitResult S4 class
#'
#' Container for the output of [tmesplit()]. Holds the decomposition matrices,
#' per-program statistics, and permutation-test p-values at the global,
#' per-group, and per-program levels.
#'
#' @slot W_shared A numeric matrix of shared-program loadings
#'   (`cell_types x k_shared`), or an empty matrix if `k_shared == 0`.
#' @slot W_specific A named list of group-specific W matrices, one per group.
#' @slot H_fractions A numeric matrix of patient-level program activities
#'   (`n_patients x k_total`), rows sum to 1 (soft assignment).
#' @slot programs A [S4Vectors::DataFrame] with one row per program
#'   (shared and specific) containing combined p-values and BH-FDR.
#' @slot p_value Numeric. Global (Level 1) p-value.
#' @slot p_value_per_group Named numeric. Per-group (Level 2) p-values.
#' @slot groups Character. Group labels in canonical order.
#' @slot k_shared Integer. Number of shared programs.
#' @slot k_specific Named integer. Number of specific programs per group.
#' @slot call Original `tmesplit()` call, for `show()`.
#' @slot metadata List of extra diagnostics (permutation settings, timing,
#'   combiner choices, etc.).
#'
#' @return An S4 object of class `TMESplitResult`. Construct an empty
#'   instance with `new("TMESplitResult")`; populated instances are
#'   returned by [tmesplit()].
#' @examples
#' res <- new("TMESplitResult")
#' isVirtualClass("TMESplitResult")
#' slotNames(res)
#' @aliases TMESplitResult
#' @export
setClass(
    "TMESplitResult",
    representation(
        W_shared = "matrix",
        W_specific = "list",
        H_fractions = "matrix",
        programs = "DataFrame",
        p_value = "numeric",
        p_value_per_group = "numeric",
        groups = "character",
        k_shared = "integer",
        k_specific = "integer",
        call = "ANY",
        metadata = "list"
    ),
    prototype(
        W_shared = matrix(numeric(0), 0, 0),
        W_specific = list(),
        H_fractions = matrix(numeric(0), 0, 0),
        p_value = NA_real_,
        p_value_per_group = numeric(0),
        groups = character(0),
        k_shared = 0L,
        k_specific = integer(0),
        metadata = list()
    )
)
