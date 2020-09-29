\name{pnormLU}
\title{Bounds for 1-Phi(.) -- Mill's Ratio related Bounds for pnorm()}
\alias{pnormL_LD10}
\alias{pnormU_S53}
\description{
  Bounds for \eqn{1 - \Phi(x)}, i.e., \code{\link{pnorm}(x, *,
    lower.tail=FALSE)}, typically related to Mill's Ratio.
}
\usage{
pnormL_LD10(x, lower.tail = FALSE, log.p = FALSE)
pnormU_S53 (x, lower.tail = FALSE, log.p = FALSE)
}
\arguments{
  \item{x}{positive (at least non-negative) numeric vector.}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{pnorm}()}.}
}
%% \details{
%% }
\value{
  a numeric vector like \code{x}
}
\references{
  Lutz Duembgen (2010)
  \emph{Bounding Standard Gaussian Tail Probabilities};
  arXiv preprint \code{1012.2063},
  \url{https://arxiv.org/abs/1012.2063}
}
\author{Martin Maechler}
%% \note{
%% }
\seealso{
 \code{\link{pnorm}}.
}
\examples{
x <- seq(1/64, 10, by=1/64)
px <- cbind(
    pnorm      (x, lower.tail=FALSE, log.p=TRUE)
  , pnormL_LD10(x, lower.tail=FALSE, log.p=TRUE)
  , pnormU_S53 (x, lower.tail=FALSE, log.p=TRUE))
matplot(x, px, type="l") # all on top of each other

matplot(x, (D <- px[,2:3] - px[,1]), type="l") # the differences
abline(h=0, lty=3, col=adjustcolor(1, 1/2))

matplot(x[x>4], D[x>4,], type="l") # the differences
abline(h=0, lty=3, col=adjustcolor(1, 1/2))

### go out !
x <- 2^seq(1, 30, by = 1/64)
px <- cbind(
    pnorm      (x, lower.tail=FALSE, log.p=TRUE)
  , pnormL_LD10(x, lower.tail=FALSE, log.p=TRUE)
  , pnormU_S53 (x, lower.tail=FALSE, log.p=TRUE))

## Shows numeric problems after ~ 10^8
matplot(x, (D <- px[,2:3] - px[,1]), type="l", log = "x")
abline(h=0, lty=3, col=adjustcolor(1, 1/2))

## Zoom into "good region":
xg <- 10 < x & x < 1e5
matplot(x[xg], D[xg,], type="l", log = "x")
abline(h=0, lty=3, col=adjustcolor(1, 1/2))

## Ok, quite revealing (?) already
xg <- 100 < x & x <= 1e4
matplot(x[xg], D[xg,], type="l", log = "x")
abline(h=0, lty=3, col=adjustcolor(1, 1/2))
}
\keyword{distribution}
