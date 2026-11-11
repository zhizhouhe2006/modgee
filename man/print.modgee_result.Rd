\name{print.modgee_result}
\alias{print.modgee_result}
\title{Print a MODGEE result}
\usage{\method{print}{modgee_result}(x, ...)}
\arguments{
\item{x}{A result returned by \code{fit_modgee()}.}
\item{...}{Additional arguments passed to \code{print()}.}
}
\description{Display the overall summary, group summaries and k comparison when present. The k comparison includes an AIC interpretation and, where supported, its direction. All other result components remain available.}
\value{The complete object, invisibly.}
