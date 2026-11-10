#' MODGEE estimation from species detection histories
#'
#' Fit the original MODGEE model to estimate detection parameters,
#' cumulative extinction and dark extinctions.
#'
#' @details
#' Use \code{\link{fit_modgee}} to fit a prepared species-by-interval
#' binary matrix. Use \code{\link{example_M}} for a documented
#' development example.
#'
#' Time aggregation and baseline selection are performed before
#' fitting. Grouped models share a temporal extinction baseline and may
#' compare fixed-k and free-k models. Bootstrap confidence intervals are
#' not computed by these functions.
#'
#' @seealso
#' \strong{Functions:}
#' \code{\link{read_matrix}}, \code{\link{fit_modgee}}, \code{\link{plot_modgee}}
#'
#' \strong{Examples:}
#' \code{\link{example_file}}, \code{\link{example_M}}, \code{\link{example_groups}}, \code{\link{example_result}}, \code{\link{example_grouped_result}}
#' @importFrom Rcpp sourceCpp
#' @useDynLib modgee, .registration = TRUE
"_PACKAGE"
