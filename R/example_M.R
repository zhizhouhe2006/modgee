#' Prepared bird detection matrix
#'
#' A numeric 0/1 matrix prepared with read_matrix(example_file, timestep = 10, baseline_min_species = 10). Species names are row names and interval labels are column names. The preparation attribute records aggregation and coverage; the final interval contains records only through 2020.
#' @format A numeric matrix.
#' @source Maintainer-provided birds_matrix.csv and birds_groups.csv.
#'   Verify original data attribution and redistribution terms before public release.
#' @seealso
#' \strong{Functions:}
#' \code{\link{read_matrix}}, \code{\link{fit_modgee}}, \code{\link{plot_modgee}}
#'
#' \strong{Examples:}
#' \code{\link{example_file}}, \code{\link{example_groups}}, \code{\link{example_result}}, \code{\link{example_grouped_result}}
#' @section Usage Examples:
#' \preformatted{
#' dim(example_M)
#' }
"example_M"
