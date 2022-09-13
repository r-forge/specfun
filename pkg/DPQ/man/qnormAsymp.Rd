\name{qnormAsymp}
\alias{qnormAsymp}
\title{Asymptotic Approximation to Outer Tail of qnorm()}
\description{
  Implementing new asymptotic tail approximations of normal quantiles,
  i.e., the \R function \code{\link{qnorm}()}, mostly useful when
  \code{log.p=TRUE} and log-scale \code{p} is relatively large negative,
  i.e., \eqn{p \ll -1}{p << -1}.
}
\usage{% ---->>>> ../R/norm_f.R <<<<--------
qnormAsymp(p,
           lp = .DT_Clog(p, lower.tail = lower.tail, log.p = log.p),
           order, lower.tail = TRUE, log.p = FALSE)
}
\arguments{
  \item{p}{numeric vector of probabilities, possibly transformed, depending
    on \code{log.p}.  Does not need to be specified, if \code{lp} is instead.}
  \item{lp}{numeric (vector) of \code{log(1-p)} values; if not specified,
    computed from \code{p}, depending on \code{lower.tail} and \code{log.p}.}
  \item{order}{an integer in \eqn{\{0,1,\dots,5\}}{{0,1,..,5}}, specifying the
    approximation order.}
  \item{lower.tail}{logical; if true, probabilities are \eqn{P[X \le x]},
    otherwise upper tail probabilities, \eqn{P[X > x]}.}
  \item{log.p}{logical; if \code{TRUE} (as typical here!), probabilities
    \eqn{p} are given as \eqn{\log(p)} in argument \code{p}.}
}
\details{
  These \emph{asymptotic} approximations have been derived by Maechler(2022)
  via iterative plug-in to the well known asymptotic approximations of
  \eqn{Q(x) = 1 - \Phi(x)} from Abramowitz and Stegun (26.2.13), p.932,
  which are provided in our package \pkg{DPQ} as \code{\link{pnormAsymp}()}.
  %%
  They will be used in R >= 4.3.0's \code{qnorm()} to provide very accurate
  quantiles in the extreme tails.
}
\value{
  a numeric vector like \code{p} or \code{lp} if that was specified instead.

  The simplemost (for extreme tails) is \code{order = 0}, where the
  asymptotic approximation is simply \eqn{\sqrt{-2s}}{sqrt(-2s)} and
  \eqn{s} is \code{-lp}.
}
\references{
  Martin Maechler (2022). Asymptotic Tail Formulas For Gaussian Quantiles;
  not yet finished.
}
\author{Martin Maechler}
\seealso{
  The upper tail approximations in Abramowitz & Stegun, in \pkg{DPQ}
  available as \code{qnormUappr()} and \code{\link{qnormUappr6}()},
  are less accurate than our \code{order >= 1} formulas in the tails.
}
\examples{%% MM: see also ./qnormUappr.Rd
pp <- c(.001, .005, .01, .05, (1:9)/10, .95, .99, .995, .999)
z_p <- qnorm(pp)
(R <- cbind(pp, z_p, qA = qnormAppr(pp)
          , qUA = qnormUappr (pp, lower.tail=TRUE)
          , qA6 = qnormUappr6(pp, lower.tail=TRUE)
          , qAsy1 = qnormAsymp(pp, lower.tail=TRUE, order=1)
          , qAsy3 = qnormAsymp(pp, lower.tail=TRUE, order=3)
          , qAsy5 = qnormAsymp(pp, lower.tail=TRUE, order=5)
          ))
## Errors, absolute and relative:
cbind(pp, (relE <- cbind(
               errA  = z_p - R[,"qA" ],
               errUA = z_p - R[,"qUA"],
               rE.A  = 1 - R[,"qA" ]/z_p,
               rE.UA = 1 - R[,"qUA"]/z_p,
               rE.A6 = 1 - R[,"qA6"]/z_p
              , rE.As1= 1 - R[,"qAsy1"]/z_p
              , rE.As3= 1 - R[,"qAsy3"]/z_p
              , rE.As5= 1 - R[,"qAsy5"]/z_p
      )))

lp <- -c(head(c(outer(c(5,2,1), 10^(4:1))), -2), 20:10, seq(9.75, 0, by = -1/8))
qnU6 <- qnormUappr6(lp=lp) # 'p' need not be specified if 'lp' is
qnAsy <- sapply(0:5, function(ord) qnormAsymp(lp=lp, lower.tail=FALSE, order=ord))
summary(warnings()) # 12 warnings ... NaNs produced
matplot(-lp, cbind(qnU6, qnAsy), type = "b", log = "x", pch=1:7)
legend("center", c("qnormUappr6()",
                paste0("qnormAsymp(*, order=",0:5,")")),
       bty="n", col=1:6, lty=1:6, pch=1:7)
matplot(-lp, cbind(qnU6, qnAsy) - qnorm(lp, lower.tail=TRUE, log.p=TRUE), pch=1:7,
        main = "absolute Error of qnorm() approximations", type = "b", log = "x")
legend("bottom", c("qnormUappr6()",
                paste0("qnormAsymp(*, order=",0:5,")")),
       bty="n", col=1:6, lty=1:6, pch=1:7)
## If you look at the numbers, in versions of R <= 4.2.x,
## qnorm() is *worse* for large -lp than the higher order approximations
}
\keyword{distribution}
