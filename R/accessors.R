#' @describeIn programs Method for [TMESplitResult].
#' @export
setMethod("programs", "TMESplitResult", function(x) x@programs)

#' @describeIn pvalue Method for [TMESplitResult]. `level = "global"` returns
#'   the Level 1 scalar, `"group"` returns the named Level 2 vector,
#'   `"program"` returns the Level 3 [S4Vectors::DataFrame].
#' @export
setMethod("pvalue", "TMESplitResult", function(x, level = "global") {
    level <- match.arg(level, c("global", "group", "program"))
    switch(level,
        global = x@p_value,
        group = x@p_value_per_group,
        program = x@programs
    )
})

#' @describeIn wShared Method for [TMESplitResult].
#' @export
setMethod("wShared", "TMESplitResult", function(x) x@W_shared)

#' @describeIn wSpecific Method for [TMESplitResult].
#' @export
setMethod("wSpecific", "TMESplitResult", function(x) x@W_specific)

#' @describeIn hFractions Method for [TMESplitResult].
#' @export
setMethod("hFractions", "TMESplitResult", function(x) x@H_fractions)
