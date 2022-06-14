\name{qtU}
\alias{qtU}
\alias{qtU1}
\title{'uniroot()'-based Computing of t-Distribution Quantiles}
\description{
  Currently, \R's own \code{\link{qt}()}  (aka \code{qnt()} in the
  non-central case) uses simple inversion of \code{\link{pt}} to compute
  quantiles in the case where \code{ncp} is specified.
  \cr
  That simple inversion (of \code{pnt()}) has seen to be deficient,
  even in cases where \code{pnt()}, i.e., \R's \code{\link{pt}(.., ncp=*)}
  does not loose accuracy.

  This \code{\link{uniroot}()}-based inversion does \emph{not} suffer from
  these deficits in some cases.
  \cr % in  ../R/t-nonc-fn.R
  \code{qtU()} is simply defined as \preformatted{qtU <- Vectorize(qtU1, c("p","df","ncp"))}
  where in \code{qtU1(p, df, ncp, *)} each of \code{(p, df, ncp)} must be of
  length one.
}
\usage{
qtU1(p, df, ncp, lower.tail = TRUE, log.p = FALSE, interval = c(-10, 10),
     tol = 1e-05, verbose = FALSE, ...)
qtU (p, df, ncp, lower.tail = TRUE, log.p = FALSE, interval = c(-10, 10),
     tol = 1e-05, verbose = FALSE, ...)
}
\arguments{
  \item{p, df, ncp}{vectors of probabilities, degrees of freedom, and
    non-centrality parameter; see \code{\link{qt}}.}
  \item{lower.tail, log.p}{logical; see \code{\link{qt}}.}
  \item{interval}{the interval in which quantiles should be searched;
    passed to \code{\link{uniroot}()};
    the current default is arbitrary and suboptimal; when \code{pt(q,*)} is
    accurate enough and hence \emph{montone} (increasing iff
    \code{lower.tail}), this interval is automatically correctly extended
    by \code{\link{uniroot}()}.}
  \item{tol}{non-negative convergence tolerance passed to \code{\link{uniroot}()}.}
  \item{verbose}{logical indicating if \emph{every} call of the objective
    function should produce a line of console output.}
  \item{\dots}{optional further arguments passed to \code{\link{uniroot}()}.}
}
\value{
  numeric vector of t quantiles, properly recycled in \code{(p, df, ncp)}.
}
\author{Martin Maechler}
\seealso{
 \code{\link{uniroot}} and \code{\link{pt}} are the simple \R level
 building blocks.  The length-1 argument version \code{qtU1()} is short and
 simple to understand.
}
\examples{
qtU1 # simple definition {with extras only for  'verbose = TRUE'}

## An example, seen to be deficient
## Stephen Berman to R-help, 13 June 2022,
## "Why does qt() return Inf with certain negative ncp values?"
q2 <- seq(-3/4, -1/4, by=1/128)
pq2 <- pt(q2, 35, ncp=-7, lower.tail=FALSE)
### ==> via qtU(), a simple uniroot() - based inversion of pt()
qpqU  <- qtU(pq2, 35, ncp=-7, lower.tail=FALSE, tol=1e-10)
stopifnot(all.equal(q2, qpqU, tol=1e-9)) # perfect!

## These two currently (2022-06-14) give Inf  whereas qtU() works fine
qt  (9e-12, df=35, ncp=-7, lower.tail=FALSE) # warnings; --> Inf
qntR(9e-12, df=35, ncp=-7, lower.tail=FALSE) #  (ditto)
## verbose = TRUE  shows all calls to pt():
qtU1(9e-12, df=35, ncp=-7, lower.tail=FALSE, verbose=TRUE)
}
