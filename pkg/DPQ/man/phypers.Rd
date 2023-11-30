\name{phypers}
\alias{phypers}
\title{The Four (4) Symmetric 'phyper()' Calls}
\description{
  Compute the four (4) symmetric \code{\link{phyper}()} calls which
  mathematically would be identical but in practice typically slightly
  differ numerically.
}
\usage{
phypers(m, n, k, q = .suppHyper(m, n, k), tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in
    \eqn{0,1,\dots, m+n}.}
  \item{q}{vector of quantiles representing the number of white balls
    drawn without replacement from an urn which contains both black and
    white balls.  By default all \dQuote{non-trivial} abscissa values
    i.e., for which the mathematical value is strictly inside \eqn{(0,1)}.
  }
  \item{tol}{a non-negative number, the \code{tolerance} for the
    \code{\link{all.equal}()} checks.}
}
\value{
  a \code{\link{list}} with components
  \item{q}{Description of 'comp1'}
  \item{phyp}{a numeric \code{\link{matrix}} of 4 columns with the 4
    different calls to \code{\link{phyper}()} which are theoretically
    equivalent because of mathematical symmetry.}
}
\references{
  Johnson et al
}
\author{Martin Maechler}
\seealso{
  \R's  \code{\link{phyper}}.  In package \CRANpkg{DPQmpfr},
  \code{\link[DPQmpfr]{phyperQ}()} uses (package \CRANpkg{gmp} based) exact
  rational arithmetic, summing up \code{\link[DPQmpfr]{dhyperQ}()}, terms
  computed by \code{\link[gmp]{chooseZ}()}, exact (long integer) arithmetic
  binomial coefficients.
}
\examples{

## The function is defined as
function(m,n,k, q = .suppHyper(m,n,k), tol = sqrt(.Machine$double.eps)) {
    N <- m+n
    pm <- cbind(ph = phyper(q,     m,  n , k), # 1 = orig.
                p2 = phyper(q,     k, N-k, m), # swap m <-> k (keep N = m+n)
                ## "lower.tail = FALSE"  <==>  1 - p..(..)
                Ip2= phyper(m-1-q, N-k, k, m, lower.tail=FALSE),
                Ip1= phyper(k-1-q, n,   m, k, lower.tail=FALSE))

    ## check that all are (approximately) the same :
    stopifnot(all.equal(pm[,1], pm[,2], tolerance=tol),
              all.equal(pm[,2], pm[,3], tolerance=tol),
              all.equal(pm[,3], pm[,4], tolerance=tol))
    list(q = q, phyp = pm)
}


str(phs <- phypers(20, 47, 31))
with(phs, cbind(q, phyp))
with(phs,
     matplot(q, phyp, type = "b"), main = "phypers(20, 47, 31)")

## differences:
with(phs, phyp[,-1] - phyp[,1])
## *relative*
relE <- with(phs, { phM <- rowMeans(phyp); 1 - phyp/phM })
print.table(cbind(q = phs$q, relE / .Machine$double.eps), zero.print = ".")
}
\keyword{arith}
