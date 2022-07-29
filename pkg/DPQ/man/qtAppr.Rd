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
  close to the Gaussian / normal distribution, using up to 4 asymptotic
  terms from Abramowitz & Stegun 26.7.5 (p.949).

}
\usage{
qtAppr (p, df, ncp, lower.tail = TRUE, log.p = FALSE, method = c("a", "b", "c"))
qtNappr(p, df,      lower.tail = TRUE, log.p=FALSE, k) % ../R/t-dist.R
}
\arguments{%--> ~/R/D/r-devel/R/src/library/stats/man/TDist.Rd : qt(p, df, .....)
  \item{p}{vector of probabilities.}
  \item{df}{degrees of freedom \eqn{> 0}, maybe non-integer.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta}; ....}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qt}()}.}
  \item{method}{a string specifying the approximation method to be used.}
  \item{k}{an integer in {0,1,2,3,4}, choosing the number of terms in \code{qtNappr()}.}
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
\examples{%% More of this >>>>> ../tests/t-nonc-tst.R
qts <- function(p, df) {
    cbind(qt = qt(p, df=df)
        , qtN0 = qtNappr(p, df=df, k=0)
        , qtN1 = qtNappr(p, df=df, k=1)
        , qtN2 = qtNappr(p, df=df, k=2)
        , qtN3 = qtNappr(p, df=df, k=3)
        , qtN4 = qtNappr(p, df=df, k=4)
          )
}
p <- (0:100)/100
ii <- 2:100 # drop p=0 & p=1  where q*(p, .) == +/- Inf

df <- 100 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
qsp1c <- qts(p, df = df)
matplot(p, qsp1c, type="l") # "all on top"
(dq <- (qsp1c[,-1] - qsp1c[,1])[ii,])
matplot(p[ii], dq, type="l", col=2:6,
        main = paste0("difference qtNappr(p,df) - qt(p,df), df=",df), xlab=quote(p))
matplot(p[ii], pmax(abs(dq), 1e-17), log="y", type="l", col=2:6,
        main = paste0("abs. difference |qtNappr(p,df) - qt(p,df)|, df=",df), xlab=quote(p))
legend("bottomright", paste0("k=",0:4), col=2:6, lty=1:5, bty="n")
matplot(p[ii], abs(dq/qsp1c[ii,"qt"]), log="y", type="l", col=2:6,
        main = sprintf("rel.error qtNappr(p, df=\%g, k=*)",df), xlab=quote(p))
legend("left", paste0("k=",0:4), col=2:6, lty=1:5, bty="n")

df <- 2000 # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
qsp1c <- qts(p, df=df)
(dq <- (qsp1c[,-1] - qsp1c[,1])[ii,])
matplot(p[ii], dq, type="l", col=2:6,
        main = paste0("difference qtNappr(p,df) - qt(p,df), df=",df), xlab=quote(p))
legend("top", paste0("k=",0:4), col=2:6, lty=1:5)
matplot(p[ii], pmax(abs(dq), 1e-17), log="y", type="l", col=2:6,
        main = paste0("abs.diff. |qtNappr(p,df) - qt(p,df)|, df=",df), xlab=quote(p))
legend("right", paste0("k=",0:4), col=2:6, lty=1:5, bty="n")

matplot(p[ii], abs(dq/qsp1c[ii,"qt"]), log="y", type="l", col=2:6,
        main = sprintf("rel.error qtNappr(p, df=\%g, k=*)",df), xlab=quote(p))
legend("left", paste0("k=",0:4), col=2:6, lty=1:5, bty="n")
}
\keyword{distribution}
