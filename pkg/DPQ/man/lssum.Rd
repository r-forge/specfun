\name{lssum}
\alias{lssum}
\title{Compute Logarithm of a Sum with Signed Large Coefficients}
\usage{
lssum(lxabs, signs, l.off = max(lxabs), strict = TRUE)
}
\arguments{
  \item{lxabs}{n-vector of values \eqn{log(|x_1|),..,log(|x_n|)}.}
  \item{signs}{corresponding signs \eqn{sign(x_1), .., sign(x_n)}.}
  \item{l.off}{the offset to substract and re-add; ideally in the order of \code{\link{max}(.)}.}
  \item{strict}{logical indicating if the function should stop on some negative sums.}
}
\value{
  \deqn{
    log(x_1 + .. + x_n) =
    = log(sum(x)) = log(sum(sign(x)*|x|)) =
    = log(sum(sign(x)*exp(log(|x|)))) =
    = log(exp(log(x0))*sum(signs*exp(log(|x|)-log(x0)))) =
    = log(x0) + log(sum(signs* exp(log(|x|)-log(x0)))) =
    = l.off + log(sum(signs* exp(lxabs - l.off )))
  }
}
\description{
  Simple vector version of \code{copula:::lssum()} (CRAN package
  \CRANpkg{copula}). %  ~/R/Pkgs/copula/R/special-func.R
  Properly compute \eqn{\log(x_1 + \ldots + x_n)}{log(x1 + .. + xn)}
  for given \eqn{log(|x_1|),.., log(|x_n|)}{log(|x_1|),.., log(|x_n|)}
  and corresponding signs \eqn{sign(x_1),.., sign(x_n)}.  Here,
  \eqn{x_i}{xi} is of arbitrary sign.
}

\seealso{
  \code{\link{lsum}()} which computes an exponential sum in log scale
  with\emph{out} signs.
}
\author{
  Marius Hofert and Martin Maechler (for package \CRANpkg{copula}).
}
\examples{
}
\keyword{distribution}
