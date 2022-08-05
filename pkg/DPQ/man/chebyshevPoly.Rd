\name{chebyshevPoly}
\title{Chebyshev Polynomial Evaluation}
\alias{chebyshevPoly}
\alias{chebyshevEval}
\alias{chebyshev_nc}
\description{
  Provides (evaluation of) Chebyshev polynomials, given their coefficients
  vector \code{coef}.
  Specifically, the following sum is evaluated:
  \deqn{\sum_{j=0}^n c_j T_j(x)}{sum(j=0..n; c[j] T[j](x))}
  where  \eqn{c_j :=}{c[j] :=} \code{coef[j+1]},
  \eqn{n :=} \code{chebyshev_nc(coef, .)} is the maximal degree and hence
  one less than the number of terms, and \eqn{T_j()}{T[j]()} is the
  Chebyshev polynomial (of the first kind) of degree \eqn{j}.
}
\usage{% --> ../R/utils.R
chebyshevPoly(coef, nc = chebyshev_nc(coef, eta), eta = .Machine$double.eps/20)

chebyshev_nc(coef, eta = .Machine$double.eps/20)
chebyshevEval(x, coef,
              nc = chebyshev_nc(coef, eta), eta = .Machine$double.eps/20)
}
\arguments{
  \item{coef}{a numeric vector of coefficients for the Chebyshev polynomial.}
  \item{nc}{the maximal degree, i.e., one less than the number of
    polynomial terms to use; typically use the default.}
  \item{eta}{a positive number; typically keep the default.}
  \item{x}{for \code{chebyshevEval()}: numeric vector of abscissa values at
    which the polynomial should be evaluated.  Typically \code{x} values
    are inside the interval \eqn{[-1, 1]}.}
}
%% \details{

%% }
\value{
  \code{chebyshevPoly()} returns \code{\link{function}(x)} which computes
  the values of the underlying Chebyshev polynomial at \code{x}.

  \code{chebyshev_nc()} returns an \code{\link{integer}}, and
  \code{chebyshevEval(x, coef)} returns a numeric \dQuote{like} \code{x}
  with the values of the polynomial at \code{x}.
}
\references{
  \url{https://en.wikipedia.org/wiki/Chebyshev_polynomials}
}
\author{R Core team, notably Ross Ihaka; Martin Maechler provided the \R interface.}
\seealso{
  \code{\link{polyn.eval}()} from CRAN package \CRANpkg{sfsmisc}; as one
  example of many more.
}
\examples{
 curve(chebyshevPoly(c(1,0,0,0))(x), -1,1, col=2, lwd=2, ylim=c(-1,1))
 abline(h=0, lty=2)
 curve(chebyshevPoly(c(0,1,0,0))(x), col=3, lwd=2, add=TRUE)
 curve(chebyshevPoly(c(0,0,1,0))(x), col=4, lwd=2, add=TRUE)
 curve(chebyshevPoly(c(0,0,0,1))(x), col=5, lwd=2, add=TRUE)

 curve(chebyshevPoly(1)(x), -1,1, col=2, lwd=2, ylim=c(-1,1))
 abline(h=0, lty=2)
 curve(chebyshevPoly(0:1)       (x), col=3, lwd=2, add=TRUE)
 curve(chebyshevPoly(c(0,0,1))  (x), col=4, lwd=2, add=TRUE)
 curve(chebyshevPoly(c(0,0,0,1))(x), col=5, lwd=2, add=TRUE)
}
\keyword{arithmetic}
