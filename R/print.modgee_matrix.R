#' Print a prepared detection matrix
#'
#' @param x A matrix returned by \code{read_matrix()}.
#' @param ... Additional arguments passed to \code{print.default()}.
#' @return The original object, invisibly. Preparation metadata are retained.
#' @export
print.modgee_matrix <- function(x, ...) {
  display <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x),
                    dimnames = dimnames(x))
  print.default(display, ...)
  invisible(x)
}
