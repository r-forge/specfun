\name{format01prec}
\alias{format01prec}
\title{Format Numbers in [0,1] with "Precise" Result}
\description{
  Format numbers in [0,1] with \dQuote{precise} result,
  notably using \code{"1-.."} if needed.
}
\usage{
format01prec(x, digits = getOption("digits"), width = digits + 2,
             eps = 1e-06, ...,
             FUN = function(x, ...) formatC(x, flag = "-", ...))
}
\arguments{
  \item{x}{numbers in [0,1]; (still works if not)}
  \item{digits}{number of digits to use; is used as
    \code{FUN(*, digits = digits)} or
    \code{FUN(*, digits = digits - 5)} depending on \code{x} or \code{eps}.}
  \item{width}{desired width (of strings in characters), is used as
    \code{FUN(*, width = width)} or
    \code{FUN(*, width = width - 2)} depending on \code{x} or \code{eps}.}
  \item{eps}{small positive number: Use \code{'1-'} for those \code{x}
    which are in \eqn{(1-eps, 1]}.  The author has claimed in the last
    millennium that (the default) 1e-6 is \emph{optimal}.}
  \item{\dots}{optional further arguments passed to \code{FUN(x, digits,
      width, ...)}.}
  \item{FUN}{a \code{\link{function}} used for \code{\link{format}()}ing;
    must accept both a \code{digits} and \code{width} argument.}
}
\value{
  a \code{\link{character}} vector of the same length as \code{x}.
}
\author{Martin Maechler, 14 May 1997}
\seealso{
  \code{\link{formatC}}, \code{\link{format.pval}}.
}
\examples{
## Show that format01prec()  does reveal more precision :
cbind(format      (1 - 2^-(16:24)),
      format01prec(1 - 2^-(16:24)))

## a bit more variety
e <- c(2^seq(-24,0, by=2), 10^-(7:1))
ee <- sort(unique(c(e, 1-e)))
noquote(ff <- format01prec(ee))
data.frame(ee, format01prec = ff)
}
\keyword{print}
\keyword{character}

