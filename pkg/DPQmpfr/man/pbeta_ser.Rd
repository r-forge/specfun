\name{pbeta_ser}
\alias{pbeta_ser}
\title{Beta Distribution Function -- \sQuote{BPSER} Series Expansion from TOMS 708}
\description{
  Compute a version of the Beta cumulative distribution function
  (\code{\link{pbeta}()} in \R), namely using the series expansion, named
  \code{BPSER()}, from \dQuote{TOMS 708}, i.e., Didonato and Morris (1992).

  This \dQuote{pure R} function exists for didactical or documentational reasons on one hand,
  as \R's own \code{\link{pbeta}()} uses this expansion when appropriate and
  other algorithms otherwise.
  On the other hand, using high precision \code{q} and MPFR arithmetic (via
  package \CRANpkg{Rmpfr}) may allow to get highly accurate \code{pbeta()}
  values.
}
\usage{
pbeta_ser(q, shape1, shape2, log.p = FALSE, eps = 1e-15, errPb = 0, verbose = FALSE)
}
\arguments{
  \item{q, shape1, shape2}{quantiles and shape parameters of the Beta
    distribution, \code{q} typically in \eqn{[0,1]}, see \code{\link{pbeta}}.
    Here, \code{q} \emph{must be scalar}, i.e., of length one, and may inherit from
    \code{class} \code{"mpfr"}, in order to be more accurate (than with the
    double precision computations).}
  \item{log.p}{if TRUE, probabilities \code{p} are given as \code{log(p)}.}
  \item{eps}{non-negative number; \code{tol <- eps/shape1} will be used for
    convergence checks in the series computations.}
  \item{errPb}{an integer code, typically in \code{-2, -1, 0} to determine how
    warnings on convergence failures are handled.}
  \item{verbose}{logical indicating if console output about intermediate
    results should be printed.}
}
\details{
  \code{pbeta_ser()} crucially needs three auxiliary functions which we
  \dQuote{mpfr-ized} as well: \code{\link{gam1M}()},
  \code{\link{lgamma1pM}()}, and \code{\link{algdivM}}.
}
\value{
  An approximation to the Beta probability  \eqn{P[X \le q]}{P[X <= q]}
  for \eqn{ X \sim B(a,b),}{ X ~ B(a,b),} (where \eqn{a=}\code{shape1}, and \eqn{b=}\code{shape2}).
}
\references{
  Didonato, A. and Morris, A., Jr, (1992)
  Algorithm 708: Significant digit computation of the incomplete beta
  function ratios,
  \emph{ACM Transactions on Mathematical Software} \bold{18}, 360--373; % =: TOMS 708
  \doi{10.1145/131766.131776}.
}
\author{
  Didonato and Morris and R Core team;
  separate packaging by Martin Maechler.
}
\seealso{
  \code{\link{pbeta}}, \pkg{DPQmpfr}'s own \code{\link{pbetaD94}};
  even more \code{pbeta()} approximations in package \CRANpkg{DPQ}, e.g.,
  \code{\link[DPQ]{pnbetaAS310}}, or \code{\link[DPQ]{pbetaRv1}}.

  In addition, for integer shape parameters, the potentially \dQuote{fully accurate}
  finite sum base \code{\link[Rmpfr]{pbetaI}()} in package \CRANpkg{Rmpfr}.
%% TODO: validate pbeta_ser() from pbetaI()
}
\examples{
(p. <- pbeta_ser(1/2, shape1 = 2, shape2 = 3, verbose=TRUE))
(lp <- pbeta_ser(1/2, shape1 = 2, shape2 = 3, log.p = TRUE))
          all.equal(lp, log(p.), tolerance=0) # 1.48e-16
stopifnot(all.equal(lp, log(p.), tolerance = 1e-13))

## Using  Vectorize() in order to allow vector 'q' e.g. for curve():
str(pbetaSer <- Vectorize(pbeta_ser, "q"))
curve(pbetaSer(x, 1.5, 4.5)); abline(h=0:1, v=0:1, lty=2, col="gray")
curve(pbeta   (x, 1.5, 4.5), add=TRUE, col = adjustcolor(2, 1/4), lwd=3)

## now using mpfr-numbers:
half <- 1/Rmpfr::mpfr(2, 256)
(p2 <- pbeta_ser(half, shape1 = 1, shape2 = 123))

}
\keyword{distribution}
\keyword{math}
