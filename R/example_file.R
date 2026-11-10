#' Path to the installed example detection CSV
#'
#' A character scalar giving the path to the bird detection CSV included in
#' the installed package. It is resolved each time the namespace is loaded.
#' The original file contains annual records from 1819 through 2020.
#' @format A character scalar containing an existing CSV path.
#' @seealso
#' \strong{Functions:}
#' \code{\link{read_matrix}}, \code{\link{fit_modgee}}, \code{\link{plot_modgee}}
#'
#' \strong{Examples:}
#' \code{\link{example_M}}, \code{\link{example_groups}}, \code{\link{example_result}}, \code{\link{example_grouped_result}}
#' @export
#' @section Usage Examples:
#' \preformatted{
#' example_file
#' file.exists(example_file)
#' }
example_file <- NULL

.onLoad <- function(libname, pkgname) {
  assign("example_file", system.file("extdata", "birds_matrix.csv",
    package = pkgname, lib.loc = libname, mustWork = TRUE),
    envir = asNamespace(pkgname))
}
