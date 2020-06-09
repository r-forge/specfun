\name{lfastchoose}
\title{R versions of Simple Formulas for Logarithmic Binomial Coefficients}
\alias{lfastchoose}
\alias{f05lchoose}
\description{
  Provide \R versions of simple formulas for computing the logarithm of
  (the absolute value of) binomial coefficients, i.e., simpler, more direct
  formulas than what (the C level) code of \R's \code{\link{lchoose}()}
  computes.
}
\usage{
lfastchoose(n, k)
 f05lchoose(n, k)
}
\arguments{
  \item{n}{a numeric vector.}
  \item{k}{a integer valued numeric vector.}
}
\value{
  a numeric vector with the same attributes as \code{n + k}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{lchoose}}.
}
\examples{
lfastchoose # function(n, k) lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
f05lchoose  # function(n, k) lfastchoose(n = floor(n + 0.5), k = floor(k + 0.5))

## interesting cases ?

}
\keyword{math}
