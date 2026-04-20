#' TMESplit discrete color palette
#'
#' A 35-color palette used as the default for categorical annotations in
#' [plotPrograms()], [plotActivities()], [plotSignificance()],
#' [plotNetwork()], and [plotSampleHeatmap()]. Colors are ordered so that the
#' first entries are visually distinct; user code can index into the vector
#' freely (e.g. `tme_palette[1:4]`) or pass a named subset as
#' `annotation_colors` to override the defaults.
#'
#' @format A character vector of 35 hex color codes.
#'
#' @examples
#' tme_palette[1:4]
#' # assign a role manually
#' sex_colors <- c(Female = tme_palette[2], Male = tme_palette[1])
#'
#' @export
tme_palette <- c(
    "#6495EDFF", "#FF69B4FF", "#00BFB2",   "#FBC4AB", "#EF8354",
    "#BA55D3FF", "#EFEA5A",   "#A4BEF3",   "#78BC61", "#F08080FF",
    "#32CD32FF", "#BB342F",   "#4357AD",   "#9ACD32FF", "#4682B4FF",
    "#DDA0DDFF", "#FFA07AFF", "#8FBC8BFF", "#4A2040", "#40E0D0FF",
    "#F0E68CFF", "#5F9EA0FF", "#D2B48CFF", "#424874", "#EF798A",
    "#FFDAB9FF", "#87CEEBFF", "#B4A0E5",   "#5BC0BE", "#773344",
    "#E0FF4F",   "#59CD90",   "#735D78",   "#FFE381", "#DE6449"
)

#' Diverging color stops for loading / z-score heatmaps
#'
#' Three-stop diverging palette (`#154999` → `white` → `#CF0034`) used by
#' [plotPrograms()] and [plotSampleHeatmap()] for W loadings and sample
#' z-scores.
#'
#' @format A character vector of length 3.
#' @export
tme_diverging <- c("#154999", "white", "#CF0034")
