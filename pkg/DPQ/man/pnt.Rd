\name{pnt}
\title{Non-central t Probability Distribution - Algorithms and Approximations}
\alias{pntR}
\alias{pntR1}
\alias{pntP94}
\alias{pntP94.1}
\alias{pntLrg}
\alias{pntJW39}
\alias{pntJW39.0}
\alias{pntChShP94}
\alias{pntChShP94.1}
\alias{pnt3150}
\alias{pnt3150.1}
\description{
 Compute different approximations for the non-central t-Distribution
 cumulative probability distribution function.
}
\usage{% --> ~/R/D/r-devel/R/src/library/stats/man/TDist.Rd
%% source in ../R/t-nonc-fn.R
pntR      (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                                          itrmax = 1000, errmax = 1e-12, verbose = TRUE)
pntR1     (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                                          itrmax = 1000, errmax = 1e-12, verbose = TRUE)
pnt3150   (t, df, ncp, lower.tail = TRUE, log.p = FALSE, M = 1000, verbose = TRUE)
pnt3150.1 (t, df, ncp, lower.tail = TRUE, log.p = FALSE, M = 1000, verbose = TRUE)
pntP94    (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                                          itrmax = 1000, errmax = 1e-12, verbose = TRUE)
pntP94.1  (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                                          itrmax = 1000, errmax = 1e-12, verbose = TRUE)
pntLrg    (t, df, ncp, lower.tail = TRUE, log.p = FALSE)
pntJW39   (t, df, ncp, lower.tail = TRUE, log.p = FALSE)
pntJW39.0 (t, df, ncp, lower.tail = TRUE, log.p = FALSE)
pntChShP94  (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                                            itrmax = 1000, errmax = 1e-12, verbose = TRUE)
pntChShP94.1(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                                            itrmax = 1000, errmax = 1e-12, verbose = TRUE)
}
\arguments{%% from Tdist.Rd (see above)
  \item{t}{vector of quantiles (called \code{q} in \code{\link{pt}(..)}.}
  %% \item{n}{number of observations. If \code{length(n) > 1}, the length
  %%   is taken to be the number required.}
  \item{df}{degrees of freedom (\eqn{> 0}, maybe non-integer).  \code{df
      = Inf} is allowed.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta};
    currently except for \code{rt()}, only for \code{abs(ncp) <= 37.62}.
    If omitted, use the central t distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}

  \item{itrmax}{ ... }
  \item{errmax}{ ... }
  \item{verbose}{ ... }
  \item{M}{positive integer specifying the number of terms to use in the series.}
}

\details{
  \describe{

    \item{\code{.}:}{ .. }
    \item{\code{.}:}{ .. }

    \item{\code{pnt3150()}, \code{pnt3150.1()}:}{
      Simple inefficient but hopefully correct version of pntP94..()
      This is really a direct implementation of formula
      (31.50), p.532 of  Johnson, Kotz and Balakrishnan (1995)
    }
    \item{\code{pntP94()}, \code{pntP94.1()}:}{New versions of
      \code{pntR1()}, \code{pntR()}; using the  Posten (1994) algorithm.
      \code{pntP94()} is the \code{\link{Vectorize}()}d version of
      \code{pntP94.1()}.}

    \item{\code{.}:}{ .. }
    \item{\code{.}:}{ .. }
  }
}

\value{
  a number for \code{pntJKBf1()} and \code{.pntJKBch1()}.

  a numeric vector of the same length as the maximum of the lengths of
  \code{x, df, ncp} for  \code{pntJKBf()} and \code{.pntJKBch()}.
}
\references{
  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995)
  Continuous Univariate Distributions Vol~2, 2nd ed.; Wiley.
  \cr Chapter 31, Section \emph{5  Distribution Function}, p.514 ff
}
\author{Martin Maechler}
\seealso{
  \code{\link{pt}}, for \R's version of non-central t probabilities.
}
\examples{
tt <- seq(0, 10, len = 21)
ncp <- seq(0, 6, len = 31)
dt3R   <- outer(tt, ncp, pt, , df = 3)
dt3JKB <- outer(tt, ncp, pntR, df = 3)# currently verbose
stopifnot(all.equal(dt3R, dt3JKB, tolerance = 4e-15))# 64-bit Lnx: 2.78e-16
}
\keyword{distribution}
\keyword{math}
