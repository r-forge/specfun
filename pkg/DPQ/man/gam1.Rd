\name{gam1}
\alias{gam1}
\title{Compute  1/Gamma(x+1) - 1  Accurately}
\description{
  Computes \eqn{1/\Gamma(a+1) - 1} accurately in \eqn{[-0.5, 1.5]} for numeric argument \code{a};
  For \code{"mpfr"} numbers, the precision is increased intermediately such
  that \eqn{a+1} should not lose precision.

 FIXME: "Pure-R" implementation is in \file{ ~/R/Pkgs/DPQ/TODO_R_versions_gam1_etc.R }
}
\usage{
gam1(a, useDbl = is.numeric(a))#, verbose = FALSE)
}
\arguments{
  \item{a}{a numeric or numeric-alike, typically inheriting from \code{class} \code{"mpfr"}.}
  \item{useDbl}{logical determining if the double precision approximation
    in \eqn{[-0.5, 1.5]} from TOMS 708 should be used.}
  %% \item{verbose}{ ..}
}
\details{ %% MM: This is all correct --- but 100% irrelevant !!!
  %%  It would be the approximation at a = -1 which is outside [-1/2, 3/2] !!!!
  \url{https://dlmf.nist.gov/} states the well-know Taylor series for
  \deqn{\frac{1}{\Gamma(z)} = \sum_{k=1}^\infty c_k z^k}{%
    1 / Gamma(z) = sum(k=1..Inf; c_k z^k)}
  with \eqn{c_1 = 1},  \eqn{c_2 = \gamma}, (Euler's gamma, \eqn{\gamma = 0.5772...}, with
  recursion \eqn{c_k = (\gamma c_{k-1} - \zeta(2) c_{k-2} ... +(-1)^k \zeta(k-1) c_1) /(k-1)}.

  Hence,  \deqn{\frac{1}{\Gamma(z+1)}    = z+1 + \sum_{k=2}^\infty c_k (z+1)^k}
          \deqn{\frac{1}{\Gamma(z+1)} -1 = z + \gamma*(z+1)^2 + \sum_{k=3}^\infty c_k (z+1)^k}
  Consequently, for \eqn{\zeta_k := \zeta(k)},
  \eqn{c_3 = (\gamma^2 - \zeta_2)/2},
  \eqn{c_4 = \gamma^3/6 - \gamma \zeta_2/2 + \zeta_3/3}.
\preformatted{
  gam <- Const("gamma", 128)
  z <- Rmpfr::zeta(mpfr(1:7, 128))
  (c3 <- (gam^2 -z[2])/2)                       # -0.655878071520253881077019515145
  (c4 <- (gam*c3 - z[2]*c2 + z[3])/3)           # -0.04200263503409523552900393488
  (c4 <- gam*(gam^2/6 - z[2]/2) + z[3]/3)
  (c5 <- (gam*c4 - z[2]*c3 + z[3]*c2 - z[4])/4) # 0.1665386113822914895017007951
  (c5 <- (gam^4/6 - gam^2*z[2] + z[2]^2/2 + gam*z[3]*4/3 - z[4])/4)
}

}
\value{
  a numeric-alike vector like \code{a}.
}
\references{
  TOMS 708, see \code{\link{pbeta}}
}
\author{Martin Maechler building on C code of TOMS 708}
%% \note{

%% }
\seealso{
 \code{\link{gamma}}.
}
\examples{
g1 <- function(u) 1/gamma(u+1) - 1
u <- seq(-.5, 1.5, by=1/16); set.seed(1); u <- sample(u) # permuted (to check logic)
%% u <- seq(-.5, 1.5, by=1/4) ## for testing
g11   <- vapply(u, gam1, 1)
gam1. <- gam1(u)
stopifnot( all.equal(g1(u), g11) )
stopifnot( identical(g11, gam1.) )

## Comparison using Rmpfr; slightly extending the [-.5, 1.5] interval:
u <- seq(-0.525, 1.525, length.out = 2001)
gam1M <- gam1(mpfr(u, 128))
relE <- asNumeric(relErrV(gam1M, gam1(u)))

plot(relE ~ u, type="l", ylim = c(-1,1) * 2.5e-15,
     main = expression("Relative Error of " ~~ gam1(u) == frac(1, Gamma(u+1)) - 1))
grid(lty = 3); abline(v = c(-.5, 1.5), col = adjustcolor(4, 1/2), lty=2, lwd=2)

## what about the direct formula -- how bad is it really ?
relED <- asNumeric(relErrV(gam1M, g1(u)))

plot(relE ~ u, type="l", ylim = c(-1,1) * 1e-14,
     main = expression("Relative Error of " ~~ gam1(u) == frac(1, Gamma(u+1)) - 1))
lines(relED ~ u, col = adjustcolor(2, 1/2), lwd = 2)
# mtext("comparing with direct formula   1/gamma(u+1) - 1")
legend("top", c("gam1(u)", "1/gamma(u+1) - 1"), col = 1:2, lwd=1:2, bty="n")
## direct is clearly *worse* , but not catastrophical

}
\keyword{math}
