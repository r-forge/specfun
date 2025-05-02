\name{dltgammaInc}
\alias{dltgammaInc}
\title{TOMS 1006 - Fast and Accurate Generalized Incomplete Gamma Function}
\description{
  This is the \code{deltagammainc} algorithm as described in Abergel and
  Moisan (2020).

  It uses one of three different algorithms to compute \eqn{G(p,x)}, their
  conveniently \emph{scaled} generalization of the (upper and lower)
  incomplete Gamma functions,

  \deqn{G(p,x) = e^{x - p\log x} \times \gamma(p,x) \ \textrm{if} x \le p or \Gamma(p,x) \textrm{otherwise},}{%
  G(p,x) = exp(x - p log(x)) * [ gamma(p,x) if x <= p  or  Gamma(p,x) otherwise],}

  for \eqn{\gamma(p,x)}{gamma(p,x)} and  \eqn{\Gamma(p,x)}{Gamma(p,x)} the
  lower and upper incomplete gamma functions, and for all
  \eqn{x \ge 0}{x >= 0} and \eqn{p > 0}.

  \code{dltgammaInc(x,y, mu, p)} := \eqn{I_{x,y}^{\mu,p}}{I[x,y](mu,p)}
  computes their (generalized incomplete gamma) integral
  \deqn{I_{x,y}^{\mu,p} := \int_x^y s^{p-1} e^{-\mu s} \; ds.}{%
        I[x,y](mu,p) := integral[x .. y] s^(p-1) e^(-mu * s) ds.}

  Current explorations (via package \CRANpkg{Rmpfr}'s
  \code{\link[Rmpfr]{igamma}()}) show that this algorithm is \emph{less
  accurate} than \R's own \code{\link{pgamma}()} at least for small \eqn{p}.
}
\usage{
dltgammaInc(x, y, mu = 1, p)
}
\arguments{
  \item{x,y, p}{numeric vectors, typically of the same length; otherwise
    recycled to common length.}
  \item{mu}{a number different from 0; must be finite; \code{1} by default.}
}
%% \details{

%% }
\value{
  a \code{\link{list}} with three components
  \item{rho, sigma}{numeric vectors of same length (the common length of
    \code{x},\code{y},.., after recycling).
    \code{ans := rho * exp(sigma)} are the computed values of \eqn{G(p,x)}.}
  \item{method}{integer code indicating the algorithm used for the
    computation of \code{(rho, sigma)}.}
}
\references{
  Rémy Abergel and Lionel Moisan (2020)
  Algorithm 1006: Fast and accurate evaluation of a generalized incomplete gamma function,
  \emph{ACM Transactions on Mathematical Software} \bold{46}(1): 1--24.
  \doi{10.1145/3365983}
}
\author{C source from Remy Abergel and Lionel Moisan, see the reference;
  original is in \file{1006.zip} at \url{https://calgo.acm.org/}.

  Interface to \R in \pkg{DPQ}: Martin Maechler
}
\seealso{
 \pkg{Rmpfr}'s \code{\link[Rmpfr]{igamma}()}, \R's \code{\link{pgamma}()}.
}
\examples{
p <- pi
x <- y <- seq(0, 10, by = 1/8)

dltg1 <- dltgammaInc(0, y,   mu = 1, p = pi); r1 <- with(dltg1, rho * exp(sigma))  ##
dltg2 <- dltgammaInc(x, Inf, mu = 1, p = pi); r2 <- with(dltg2, rho * exp(sigma))
str(dltg1)
stopifnot(identical(c(rho = "double", sigma = "double", method = "integer"),
                    sapply(dltg1, typeof)))
summary(as.data.frame(dltg1))
pg1 <- pgamma(q = y, shape = pi, lower.tail=TRUE)
pg2 <- pgamma(q = x, shape = pi, lower.tail=FALSE)

stopifnot(all.equal(r1 / gamma(pi), pg1, tolerance = 1e-14)) # seen 5.26e-15
stopifnot(all.equal(r2 / gamma(pi), pg2, tolerance = 1e-14)) #  "   5.48e-15

if(requireNamespace("Rmpfr")) withAutoprint({
    asN  <- Rmpfr::asNumeric
    mpfr <- Rmpfr::mpfr
    iGam <- Rmpfr::igamma(a= mpfr(pi, 128),
                          x= mpfr(y,  128))
    relErrV <- sfsmisc::relErrV
    relE <- asN(relErrV(target = iGam, current = r2))
    summary(relE)
    plot(x, relE, type = "b", col=2, # not very good: at  x ~= 3, goes up to  1e-14 !
         main = "rel.Error(incGamma(pi, x)) vs  x")
    abline(h = 0, col = adjustcolor("gray", 0.75), lty = 2)
    lines(x, asN(relErrV(iGam, pg2 * gamma(pi))), col = adjustcolor(4, 1/2), lwd=2)
    ## ## ==> R's pgamma() is *much* more accurate !!  ==> how embarrassing for TOMS 1006 !?!?
    legend("topright", expression(dltgammaInc(x, Inf, mu == 1, p == pi),
                                  pgamma(x, pi, lower.tail== F) \%*\% Gamma(pi)), %
           col = c(palette()[2], adjustcolor(4, 1/2)), lwd = c(1, 2), pch = c(1, NA), bty = "n")
})

}
\keyword{math}
