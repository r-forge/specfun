\name{lgammaP11}
\alias{lgammaP11}
\title{Log Gamma(p) for Positive \code{p} by Pugh's Method (11 Terms)}
\description{
  Computes \eqn{\log(\Gamma(p))} for \eqn{p > 0},
  where \eqn{\Gamma(p) = \int_0^\infty s^{p-1} e^{-s}\; ds}{%
    Gamma(p) := Int[0, Inf]s^(p-1) exp(-s)}.

  The evaluation hapenss via Pugh's method (approximation with 11 terms),
  which is a refinement of the Lanczos method (with 6 terms).

  This implementation has been used in and for Abergel and Moisan (2020)'s
  TOMS 1006 algorithm, see the reference and our \code{\link{dltgammaInc}()}.
  Given the examples below and \code{example(dltgammaInc)}, their use of
  \dQuote{accurate} does not apply to \R's \dQuote{standard} accuracy.
}
\usage{
lgammaP11(x)
}
\arguments{
  \item{x}{a numeric vector}
}
\value{
  a numeric vector "as" \code{x} (with no attributes, currently).
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
  \code{\link{dltgammaInc}()} partly uses (the underly C code of) \code{lgammaP11()};
  \pkg{DPQ}'s \code{\link{lgamma1p}}, \R's \code{\link{lgamma}()}.
}
\examples{
(r <- lgammaP11(1:20))
stopifnot(all.equal(r, lgamma(1:20)))
summary(rE <- sfsmisc::relErrV(r, lgamma(1:20)))
cbind(x=1:20, r, lgamma(1:20))
## at x=1 and 2 it is *not* fully accurate but gives -4.44e-16  (= - 2*EPS )

if(requireNamespace("Rmpfr")) withAutoprint({
    asN  <- Rmpfr::asNumeric
    mpfr <- Rmpfr::mpfr
    x <- seq(1, 20, by = 1/8)
    lgamM <- lgamma(x = mpfr(x, 128)) # uses Rmpfr::Math(.) -> Rmpfr:::.Math.codes
    lgam  <- lgammaP11(x)
    relErrV <- sfsmisc::relErrV
    relE <- asN(relErrV(target = lgamM, current = lgam))
    summary(relE)
    summary(relE.R <- asN(relErrV(lgamM, base::lgamma(x))))
    plot(x, relE, type = "b", col=2, # not very good: at  x ~= 3, goes up to  1e-14 !
         main = "rel.Error(lgammaP11(x)) vs  x")
    abline(h = 0, col = adjustcolor("gray", 0.75), lty = 2)
    lines(x, relE.R, col = adjustcolor(4, 1/2), lwd=2)
    ## ## ==> R's lgamma() is much more accurate
    legend("topright", expression(lgammaP11(x), lgamma(x)),
           col = c(palette()[2], adjustcolor(4, 1/2)), lwd = c(1, 2), pch = c(1, NA), bty = "n")

    ## |absolute rel.Err| on log-scale:
    cEps <- 2^-52 ; lc <- adjustcolor("gray", 0.75)
    plot(x, abs(relE), type = "b", col=2, log = "y", ylim = cEps * c(2^-4, 2^6), yaxt = "n",
         main = "|rel.Error(lgammaP11(x))| vs  x -- log-scale")
    sfsmisc::eaxis(2, cex.axis = 3/4, col.axis = adjustcolor(1, 3/4))
    abline(h= cEps, col = lc, lty = 2); axis(4, at=cEps, quote(epsilon[C]), las=1, col.axis=lc)
    lines(x, pmax(2^-6*cEps, abs(relE.R)), col = adjustcolor(4, 1/2), lwd=2)
    ## ## ==> R's lgamma() is *much* more accurate
    legend("topright", expression(lgammaP11(x), lgamma(x)), bg = adjustcolor("gray", 1/4),
           col = c(palette()[2], adjustcolor(4, 1/2)), lwd = c(1, 2), pch = c(1, NA), bty = "n")
    mtext(sfsmisc::shortRversion(), cex = .75, adj = 1)
})
}
\keyword{math}
