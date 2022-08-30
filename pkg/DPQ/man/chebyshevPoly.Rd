\name{chebyshevPoly}
\title{Chebyshev Polynomial Evaluation}
\alias{chebyshevPoly}
\alias{chebyshevEval}
\alias{chebyshev_nc}
\description{
  Provides (evaluation of) Chebyshev polynomials, given their coefficients
  vector \code{coef} (using \eqn{2 c_0}, i.e., \code{2*coef[1]} as the base
  R mathlib \code{chebyshev*()} functions.
  Specifically, the following sum is evaluated:
  \deqn{\sum_{j=0}^n c_j T_j(x)}{sum(j=0..n; c[j] T[j](x))}
  where  \eqn{c_0 :=}{c[0] :=}\code{coef[1]}  and
         \eqn{c_j :=}{c[j] :=}\code{coef[j+1]} for \eqn{j \ge 1}{j >= 1}.
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
## The first 5 (base) Chebyshev polynomials:
T0 <- chebyshevPoly(2)  # !! 2, not 1
T1 <- chebyshevPoly(0:1)
T2 <- chebyshevPoly(c(0,0,1))
T3 <- chebyshevPoly(c(0,0,0,1))
T4 <- chebyshevPoly(c(0,0,0,0,1))
curve(T0(x), -1,1, col=1, lwd=2, ylim=c(-1,1))
abline(h=0, lty=2)
curve(T1(x), col=2, lwd=2, add=TRUE)
curve(T2(x), col=3, lwd=2, add=TRUE)
curve(T3(x), col=4, lwd=2, add=TRUE)
curve(T4(x), col=5, lwd=2, add=TRUE)

(Tv <- vapply(c(T0=T0, T1=T1, T2=T2, T3=T3, T4=T4),
              function(Tp) Tp(-1:1), numeric(3)))
x <- seq(-1,1, by = 1/64)
stopifnot(exprs = {
   all.equal(chebyshevPoly(1:5)(x),
             0.5*T0(x) + 2*T1(x) + 3*T2(x) + 4*T3(x) + 5*T4(x))
   all.equal(unname(Tv), rbind(c(1,-1), c(1:-1,0:1), rep(1,5)))# warning on rbind()
})

}
\keyword{arithmetic}
