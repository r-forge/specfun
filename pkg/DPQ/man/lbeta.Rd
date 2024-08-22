\name{lbeta}
\title{(Log) Beta and Ratio of Gammas Approximations}
\alias{lbetaM}
\alias{lbeta_asy}
\alias{lbetaMM}
\alias{lbetaI}
\alias{betaI}
\alias{logQab_asy}
\alias{Qab_terms}
\description{
  Compute \code{log(\link{beta}(a,b))} in a simple (fast) or asymptotic
  way.  The asymptotic case is based on the asymptotic \eqn{\Gamma}
  (\code{\link{gamma}}) ratios, provided in \code{Qab_terms()} and
  \code{logQab_asy()}.

  \code{lbeta_asy(a,b, ..)} is simply \code{lgamma(a) - logQab_asy(a, b, ..)}.
}
\usage{%% in ../R/beta-fns.R <<<<<<<
lbetaM   (a, b, k.max = 5, give.all = FALSE)
lbeta_asy(a, b, k.max = 5, give.all = FALSE)
lbetaMM  (a, b, cutAsy = 1e-2, verbose = FALSE)

 betaI(a, n)
lbetaI(a, n)

logQab_asy(a, b, k.max = 5, give.all = FALSE)
Qab_terms(a, k)
}
\arguments{
  \item{a, b, n}{the Beta parameters, see \code{\link{beta}}; \code{n} must
  be a positive integer and \dQuote{small}.}
  \item{k.max, k}{for \code{lbeta*()} and \code{logQab_asy()}: the number
    of terms to be used in the series expansion of \code{Qab_terms()},
    currently must be in \eqn{{0, 1, .., 5}}.}
  \item{give.all}{\code{\link{logical}} indicating if all terms should be
    returned (as columns of a matrix) or just the result.}
  \item{cutAsy}{cutoff value from where to switch to asymptotic formula.}%% FIXME details
  \item{verbose}{logical (or integer) indicating if and how much monitoring
    information should be printed to the console.}
}
\details{
  All \code{lbeta*()} functions compute \code{log(\link{beta}(a,b))}.

  We use \eqn{Qab = Qab(a,b)} for
  \deqn{Q_{a,b} := \frac{\Gamma(a + b)}{\Gamma(b)},}{\Gamma(a + b) / \Gamma(b),}
  which is numerically challenging when \eqn{b} becomes large compared to
  \code{a}, or \eqn{a \ll b}{a << b}.

  With the beta function
  \deqn{B(a,b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)} = \frac{\Gamma(a)}{Qab},}{%
        B(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b) = \Gamma(a) / Qab,}
  and hence
  \deqn{\log B(a,b) = \log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b) = \log\Gamma(a) - \log Qab,}{%
         log B(a,b) =  log\Gamma(a) +  log\Gamma(b) -  log\Gamma(a+b) =  log\Gamma(a) -  log Qab,}
  or in \R, \code{lbeta(a,b) := lgamma(a) - logQab(a,b)}.

  Indeed, typically everything has to be computed in log scale, as both \eqn{\Gamma(b)}
  and \eqn{\Gamma(a+b)} would overflow numerically for large \eqn{b}.
  Consequently, we use \code{logQab*()}, and for the large \eqn{b} case
  \code{logQab_asy()} specifically,
  \deqn{\code{logQab(a,b)} := \log( Qab(a,b) ).}{\code{logQab(a,b)} := log(Qab(a,b)).}
  The 5 polynomial terms in \code{Qab_terms()} have been derived by the
  author in 1997, but not published, about getting asymptotic formula for
  \eqn{\Gamma} ratios, related to but \emph{different} than formula
  (6.1.47) in Abramowitz and Stegun. % see ../R/beta-fns.R

  We also have a vignette % ../vignettes/comp-beta.Rnw
  about this, but really the problem has been adressed pragmatically
  by the authors of TOMS 708, see the \sQuote{References} in
  \code{\link{pbeta}},% ~/R/D/r-devel/R/src/nmath/toms708.c
  by their routine \code{\link{algdiv}()} which also is available in our
  package \pkg{DPQ}, \eqn{\code{algdiv}(a,b) = - \code{logQab}(a,b)}.
  Note that this is related to computing \code{\link{qbeta}()} in boundary
  cases.  See also \code{\link{algdiv}()} \sQuote{Details}.
}

\value{
  a fast or simple (approximate) computation of \code{\link{lbeta}(a,b)}.
}
\references{
  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain;
  Formula (6.1.47), p.257
}
\author{Martin Maechler}
%% \note{
%% }
\seealso{
  \R's \code{\link{beta}} function; \code{\link{algdiv}()}.
}
\examples{
(r  <- logQab_asy(1, 50))
(rF <- logQab_asy(1, 50, give.all=TRUE))
r == rF # all TRUE:  here, even the first approx. is good!
(r2  <- logQab_asy(5/4, 50))
(r2F <- logQab_asy(5/4, 50, give.all=TRUE))
r2 == r2F # TRUE only first entry "5"
(r2F.3 <- logQab_asy(5/4, 50, k=3, give.all=TRUE))

## Check relation to Beta(), Gamma() functions:
a <- 1.1 * 2^(-6:4)
b <- 1001.5
rDlgg <- lgamma(a+b) - lgamma(b) # suffers from cancellation for small 'a'
rDlgb <- lgamma(a) - lbeta(a, b) #    (ditto)
ralgd <- - algdiv(a,b)
rQasy <- logQab_asy(a, b)
cbind(a, rDlgg, rDlgb, ralgd, rQasy)
all.equal(rDlgg, rDlgb, tolerance = 0) # 3.0e-14
all.equal(rDlgb, ralgd, tolerance = 0) # 1.2e-16
all.equal(ralgd, rQasy, tolerance = 0) # 4.1e-10
all.equal(rQasy, rDlgg, tolerance = 0) # 3.5e-10

stopifnot(exprs = {
    all.equal(rDlgg, rDlgb, tolerance = 1e-12) # 3e-14 {from cancellations!}
    all.equal(rDlgb, ralgd, tolerance = 1e-13) # 1e-16
    all.equal(ralgd, rQasy, tolerance = 2e-9) # 4.1e-10
    all.equal(rQasy, rDlgg, tolerance = 2e-9) # 3.5e-10
    all.equal(lgamma(a)-lbeta(a, 2*b), logQab_asy(a, 2*b), tolerance =1e-10)# 1.4e-11
    all.equal(lgamma(a)-lbeta(a, b/2), logQab_asy(a, b/2), tolerance = 1e-7)# 1.2e-8
})
if(requireNamespace("Rmpfr")) withAutoprint({
  aM <- Rmpfr::mpfr(a, 512)
  bM <- Rmpfr::mpfr(b, 512)
  rT <- lgamma(aM+bM) - lgamma(bM) # "True" i.e. accurate values
  relE <- Rmpfr::asNumeric(sfsmisc::relErrV(rT, cbind(rDlgg, rDlgb, ralgd, rQasy)))
  cbind(a, signif(relE,4))
  ##          a      rDlgg      rDlgb      ralgd      rQasy
  ##  0.0171875  4.802e-12  3.921e-16  4.145e-17 -4.260e-16
  ##  0.0343750  1.658e-12  1.509e-15 -1.011e-17  1.068e-16
  ##  0.0687500 -2.555e-13  6.853e-16 -1.596e-17 -1.328e-16
  ##  0.1375000  1.916e-12 -7.782e-17  3.905e-17 -7.782e-17
  ##  0.2750000  1.246e-14  7.001e-17  7.001e-17 -4.686e-17
  ##  0.5500000 -2.313e-13  5.647e-17  5.647e-17 -6.040e-17
  ##  1.1000000 -9.140e-14 -1.298e-16 -1.297e-17 -1.297e-17
  ##  2.2000000  9.912e-14  2.420e-17  2.420e-17 -9.265e-17
  ##  4.4000000  1.888e-14  6.810e-17 -4.873e-17 -4.873e-17
  ##  8.8000000 -7.491e-15  1.004e-16 -1.638e-17 -4.118e-13
  ## 17.6000000  2.222e-15  1.207e-16  3.974e-18 -6.972e-10

## ==>  logQab_asy() is very good _here_ as long as  a << b
})

}%
\keyword{math}

