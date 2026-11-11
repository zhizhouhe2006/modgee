\name{print.modgee_matrix}
\alias{print.modgee_matrix}
\title{Print a prepared detection matrix}
\usage{\method{print}{modgee_matrix}(x, ...)}
\arguments{
\item{x}{A matrix returned by \code{read_matrix()}.}
\item{...}{Additional arguments passed to \code{print.default()}.}
}
\description{Display matrix entries and row/column names without preparation metadata.}
\value{The original object, invisibly. Preparation metadata are retained.}
