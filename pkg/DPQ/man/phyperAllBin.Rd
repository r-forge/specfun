\name{phyperAllBin}
\title{Compute Hypergeometric Probabilities via Binomial Approximations}
\alias{phyperAllBin}
\alias{phyperAllBinM}
\alias{.suppHyper}
\description{
  \itemize{
    \item \code{phyperAllBinM()} computes all four Molenaar binomial approximations
    to the hypergeometric cumulative distribution function \code{\link{phyper}()}.
    \item \code{phyperAllBin()}  computes Molenaar's four, plus the other four
    \code{\link{phyperBin.1}()}, \code{*.2}, \code{*.3}, and \code{*.4}.
  }
}
\usage{
phyperAllBin (m, n, k, q = .suppHyper(m, n, k), lower.tail = TRUE, log.p = FALSE)
phyperAllBinM(m, n, k, q = .suppHyper(m, n, k), lower.tail = TRUE, log.p = FALSE)
.suppHyper(m, n, k)
}
\arguments{
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in \eqn{0,1,\dots, m+n}.}
  \item{q}{vector of quantiles representing the number of white balls
    drawn without replacement from an urn which contains both black and
    white balls.  The default, \code{.suppHyper(m, n, k)} provides the full
    (finite) support.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
}
\value{
  the \code{phyperAllBin*()} functions return
  a numeric \code{\link{matrix}}, with each column a different
  approximation to \code{\link{phyper}(m,n,k,q, lower.tail, log.p)}.

  Note that the columns of \code{phyperAllBinM()} are a \emph{subset} of
  those from \code{phyperAllBin()}.
}
\references{See those in \code{\link{phyperBinMolenaar}}.}
\author{Martin Maechler}
\seealso{
 \code{\link{phyperBin.1}} etc, and
 \code{\link{phyperBinMolenaar}}.

 \code{\link{phyper}}
}
\examples{
.suppHyper # very simple:
stopifnot(identical(.suppHyper, ignore.environment = TRUE,
         function (m, n, k) max(0, k-n):min(k, m)))

phBall <- phyperAllBin (5,15, 7)
phBalM <- phyperAllBinM(5,15, 7)
stopifnot(identical(
            phBall[, colnames(phBalM)] ,
            phBalM)
         , .suppHyper(5, 15, 7) == 0:5
)

round(phBall, 4)
## relative Error: number of correct digits =
cbind(q = 0:5, round(-log10(abs(1 - phBall / phyper(0:5, 5,15,7))),  digits=2))
}
\keyword{ distribution }
