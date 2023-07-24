\name{phyperR}
\alias{phyperR}
\title{\R-only version of \R's original phyper() algorithm}
\description{
  An \R version of the first \code{phyper()} algorithm in \R, which was
  used up to svn rev \code{30227} on 2004-07-09.
}
\usage{% ../R/hyper-dist.R (at end)
phyperR(q, m, n, k, lower.tail=TRUE, log.p=FALSE)
}
\arguments{
  \item{q}{vector of quantiles representing the number of white balls
    drawn without replacement from an urn which contains both black and
    white balls.}
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in \eqn{0,1,\dots, m+n}.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
%% \details{

%% }
\note{
  The original argument list in \code{C} was \code{(x, NR, NB, n)} where
  there were \emph{red} and \emph{black} balls in the urn.

  Note that we have \emph{vectorized} a translation to \R of the original C
  code.
}
\value{
  a numeric vector similar to \code{phyper(q, m, n, k)}.
}
%% \references{
%% }
\author{Martin Maechler}
\seealso{
 \code{\link{phyper}} and our \code{\link{phyperR2}()} for the pure \R
  version of the newer (Welinder) \code{phyper()} algorithm
}
\examples{
m <- 9:12; n <- 7:10; k <- 10
x <- 0:(k+1) # length 12
## confirmation that recycling + lower.tail, log.p now work:
for(lg in c(FALSE,TRUE))
  for(lt in c(FALSE, TRUE)) {
    cat("(lower.tail = ", lt, " -- log = ", lg,"):\n", sep="")
    withAutoprint({
      (rr <-
           cbind(x, m, n, k, # recycling (to 12 rows)
                 ph  = phyper (x, m, n, k, lower.tail=lt, log.p=lg),
                 phR = phyperR(x, m, n, k, lower.tail=lt, log.p=lg)))
      all.equal(rr[,"ph"], rr[,"phR"], tol = 0)
      ## saw   4.706e-15 1.742e-15 7.002e-12 1.086e-15  [x86_64 Lnx]
      stopifnot(all.equal(rr[,"ph"], rr[,"phR"],
                          tol = if(lg && !lt) 2e-11 else 2e-14))
    })
  }

}
\keyword{ distribution }
