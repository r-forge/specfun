\name{lbeta}
\title{(Log) Beta Approximations}
\alias{lbetaM}
\alias{lbeta_asy}
\alias{lbetaMM}
\alias{lbetaI}
\alias{betaI}
\alias{logQab_asy}
\alias{Qab_terms}
\description{
  Compute \code{log(\link{beta}(a,b))} in a simple (fast) or asymptotic way.
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
  \item{k.max}{ .. }
  \item{give.all}{\code{\link{logical}} .. }
  \item{cutAsy}{cutoff value from where to switch to asymptotic formula.}%% FIXME details
  \item{verbose}{logical (or integer) indicating if and how much monitoring
    information should be printed to the console.}
  \item{k}{the number of terms in the series expansion of
    \code{Qab_terms()}, currently must be in \eqn{{0, 1, .., 5}}.}
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
        log  B(a,b) =  log\Gamma(a) + log\Gamma(b)  -  log\Gamma(a+b) = log\Gamma(a) -  log Qab,}
  or in \R, \code{lBeta(a,b) := lgamma(a) - logQab(a,b)}.

  Indeed, typically everything has to be computed in log scale, as both \eqn{\Gamma(b)}
  and \eqn{\Gamma(a+b)} would overflow numerically for large \eqn{b}.
  Consequently, we use \code{logQab*()}, and for the large \eqn{b} case
  \code{logQab_asy()} specifically,
  \deqn{\code{logQab(a,b)} := \log( Qab(a,b) ).}{\code{logQab(a,b)} := log(Qab(a,b)).}

  Note this is related to trying to get asymptotic formula for
  \eqn{\Gamma} ratios, notably formula (6.1.47) in Abramowitz and Stegun.

  Note how this is related to computing \code{\link{qbeta}()} in boundary
  cases, and see \code{\link{algdiv}()} \sQuote{Details} about this.


  We also have a vignette % ../vignettes/comp-beta.Rnw
  about this, but really the problem has been adressed pragmatically
  by the authors of TOMS 708, see the \sQuote{References} in
  \code{\link{pbeta}},% ~/R/D/r-devel/R/src/nmath/toms708.c
  by their routine \code{\link{algdiv}()} which also is available in our
  package \pkg{DPQ}.
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
  ## TODO
}
\keyword{math}

