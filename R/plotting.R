#' @describeIn plotPrograms Method for [TMESplitResult].
#' @export
setMethod("plotPrograms", "TMESplitResult", function(x, ...) {
    .plot_not_implemented("plotPrograms")
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
