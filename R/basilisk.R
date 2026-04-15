#' Basilisk environment for TMESplit
#'
#' Sandboxed conda env containing `tmesplitpy` (installed from its GitHub
#' `main` branch during Phase 3.2 development; will switch to a pinned
#' PyPI version for Bioconductor submission). Users never interact with
#' this directly — [tmesplit()] and friends call it via
#' [basilisk::basiliskRun].
#'
#' @keywords internal
#' @noRd
.tmesplit_env <- basilisk::BasiliskEnvironment(
    envname = "tmesplit_env",
    pkgname = "TMESplit",
    packages = c(
        "python=3.11",
        "numpy=1.26",
        "scipy=1.13",
        "scikit-learn=1.5",
        "pandas=2.2",
        "anndata=0.10",
        "numba=0.60",
        "joblib=1.4"
    ),
    pip = c(
        "tmesplitpy @ git+https://github.com/petrenkokate/TMESplitpy.git@main"
    )
)

#' Run a function inside the TMESplit basilisk env
#'
#' Thin wrapper around [basilisk::basiliskRun] that pins the environment
#' so callers don't have to.
#'
#' @param fun Function evaluated inside the env. Receives no positional
#'   args; use a closure to capture data.
#' @param ... Passed to `fun`.
#' @keywords internal
#' @noRd
.basilisk_run <- function(fun, ...) {
    proc <- basilisk::basiliskStart(.tmesplit_env)
    on.exit(basilisk::basiliskStop(proc), add = TRUE)
    basilisk::basiliskRun(proc, fun, ...)
}
