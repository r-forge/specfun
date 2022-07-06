\name{qtAppr}
\alias{qtAppr}
\alias{qtNappr}
\title{Compute Approximate Quantiles of the (Non-Central) t-Distribution}
\description{
  Compute quantiles (inverse distribution values) for the non-central t distribution.
  using Johnson,Kotz,.. p.521, formula (31.26 a) (31.26 b) & (31.26 c)

  Note that \code{\link{qt}(.., ncp=*)} did not exist yet in 1999, when MM
  implemented \code{qtAppr()}.

  \code{qtNappr()} approximates t-quantiles for large \code{df}, i.e., when
  close to the Gaussian / normal distribution, using Abramowitz \& Stegun 26.7.5 (p.949).

}
\usage{
qtAppr (p, df, ncp, lower.tail = TRUE, log.p = FALSE, method = c("a", "b", "c"))
qtNappr(p, df,      lower.tail=TRUE, log.p=FALSE, k = 2) % ../R/t-dist.R
}
\arguments{%--> ~/R/D/r-devel/R/src/library/stats/man/TDist.Rd : qt(p, df, ncp, lower.tail = TRUE, log.p = FALSE)
  \item{p}{vector of probabilities.}
  \item{df}{degrees of freedom \eqn{> 0}, maybe non-integer.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta}; ....}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qt}()}.}
  \item{method}{a string specifying the approximation method to be used.}
  \item{k}{an integer in {0,1,2}, choosing the number of terms in \code{qtNappr()}.}
}
%% \details{
%% }
\value{
  numeric vector of length \code{length(p + df + ncp)} with approximate t-quantiles.
}
\references{
  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995)
  Continuous Univariate Distributions Vol~2, 2nd ed.; Wiley.
  \cr Chapter 31, Section \emph{6  Approximation}, p.519 ff

  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \cr Formula (26.7.5), p.949
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain.
}
%\author{Martin Maechler}
\seealso{
  Our \code{\link{qtU}()}; several non-central density and probability
  approximations in \code{\link{dntJKBf}}, and e.g., \code{\link{pntR}}.
  Further, \R's \code{\link{qt}}.
}
\examples{%% TODO move more of this to bottom of  ../tests/t-nonc-tst.R
  qts <- function(p, df, ...) {
      cbind(qt = qt(p, df=df)
          , qtN0 = qtNappr(p, df=df,k=0)
          , qtN1 = qtNappr(p, df=df,k=1)
          , qtN2 = qtNappr(p, df=df,k=2)
            )
  }
  p <- (0:100)/100
  ##  df = 100
  qsp1c <- qts(p, df = 100)
  matplot(p, qsp1c, type="l") # "all on top"
  (dq <- (qsp1c[,-1] - qsp1c[,1])[2:100,])
  matplot(p[2:100], dq, type="l", col=2:4,
          main = "difference qtNappr(p,df) - qt(p,df), df=100", xlab=quote(p))
  legend("top", paste0("k=",0:2), col=2:4, lty=1:3)
  ##  df = 2000
  qsp1c <- qts(p, df = 2000)
  matplot(p, qsp1c, type="l") # "all on top"
  (dq <- (qsp1c[,-1] - qsp1c[,1])[2:100,])
  matplot(p[2:100], dq, type="l", col=2:4,
          main = "difference qtNappr(p,df) - qt(p,df), df=2000", xlab=quote(p))
  legend("top", paste0("k=",0:2), col=2:4, lty=1:3)
  ## drop the k=0 one
  matplot(p[2:100], dq[,-1], type="l", col=3:4,
          main = "difference qtNappr(p,df) - qt(p,df), df=2000", xlab=quote(p))
  legend("top", paste0("k=",1:2), col=3:4, lty=1:3)
  ## k=2 for df=2000:
  plot(p[2:100], dq[,3], type="l", main="absolute difference .. k=2, df=2000")
  plot(p[2:100], dq[,3]/qsp1c[2:100,1], type="l", main="relative difference .. k=2, df=2000")
  ## qtAppr is really worse (than k=2 (?) --- FIXME: look closer!
  qtA <- qtAppr(p, df=2000, ncp=0)
  lines(p, qtA/qsp1c[,1] - 1, col=2, lwd=2, lty=2)
  cbind(p, cbind(qsp1c, qtA)) ## qtA  quite close to k=2 ???
  cbind(p, cbind(qsp1c[,-1], qtA) - qsp1c[,1]) ## qtA  between k=1 and k=2 ...
}
\keyword{distribution}
