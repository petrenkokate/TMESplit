#' Show method for TMESplitResult
#' @param object A [TMESplitResult].
#' @return Invisibly returns `NULL`; called for its side effect (printing).
#' @examples
#' show(new("TMESplitResult"))
#' @export
setMethod("show", "TMESplitResult", function(object) {
    cat("TMESplitResult\n")
    cat("  groups:     ",
        paste(object@groups, collapse = ", "), "\n", sep = "")
    cat("  k_shared:   ", object@k_shared, "\n", sep = "")
    if (length(object@k_specific)) {
        ks <- paste(names(object@k_specific),
                    unname(object@k_specific), sep = "=",
                    collapse = ", ")
        cat("  k_specific: ", ks, "\n", sep = "")
    }
    if (length(object@p_value) && !is.na(object@p_value)) {
        cat("  p (global):      ",
            format.pval(object@p_value, digits = 3),
            "\n", sep = "")
    }
    if (length(object@p_value_per_group)) {
        pg <- paste(names(object@p_value_per_group),
                    format.pval(object@p_value_per_group, digits = 3),
                    sep = "=", collapse = ", ")
        cat("  p (per group):   ", pg, "\n", sep = "")
    }
    if (length(object@programs)) {
        cat("  programs:   ", nrow(object@programs),
            " (use programs(x) to inspect)\n", sep = "")
    }
    invisible(NULL)
})
