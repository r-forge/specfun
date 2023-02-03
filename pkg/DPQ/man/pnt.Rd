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
           use.pnorm = (df > 4e5 ||
                        ncp^2 > 2*log(2)*(-.Machine$double.min.exp)),
                                          itrmax = 1000, errmax = 1e-12, verbose = TRUE)
pntR1     (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
           use.pnorm = (df > 4e5 ||
                        ncp^2 > 2*log(2)*(-.Machine$double.min.exp)),
                                          itrmax = 1000, errmax = 1e-12, verbose = TRUE)

pntP94    (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                                          itrmax = 1000, errmax = 1e-12, verbose = TRUE)
pntP94.1  (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
                                          itrmax = 1000, errmax = 1e-12, verbose = TRUE)

pnt3150   (t, df, ncp, lower.tail = TRUE, log.p = FALSE, M = 1000, verbose = TRUE)
pnt3150.1 (t, df, ncp, lower.tail = TRUE, log.p = FALSE, M = 1000, verbose = TRUE)

pntLrg    (t, df, ncp, lower.tail = TRUE, log.p = FALSE)

pntJW39   (t, df, ncp, lower.tail = TRUE, log.p = FALSE)
pntJW39.0 (t, df, ncp, lower.tail = TRUE, log.p = FALSE)
% TODO: *IMPLEMENT* (and then document in \details{}) !
%% pntChShP94  (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
%%                                             itrmax = 1000, errmax = 1e-12, verbose = TRUE)
%% pntChShP94.1(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
%%                                             itrmax = 1000, errmax = 1e-12, verbose = TRUE)
}
\arguments{%% from Tdist.Rd (see above)
  \item{t}{vector of quantiles (called \code{q} in \code{\link{pt}(..)}).}
  \item{df}{degrees of freedom (\eqn{> 0}, maybe non-integer).  \code{df
      = Inf} is allowed.}
  \item{ncp}{non-centrality parameter \eqn{\delta \ge 0};
    If omitted, use the central t distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}

  \item{use.pnorm}{\code{\link{logical}} indicating if the
    \code{\link{pnorm}()} approximation of Abramowitz and Stegun (26.7.10)
    should be used, which is available as \code{pntLrg()}.

    The default corresponds to \R \code{\link{pt}()}'s own behaviour (which
    is most probably suboptimal).}

  \item{itrmax}{number of iterations / terms.}% .. say more ?
  \item{errmax}{convergence bound for the iterations.}%% could say more
  \item{verbose}{\code{\link{logical}} or integer determining the amount of
    diagnostic print out to the console.}
  \item{M}{positive integer specifying the number of terms to use in the series.}
}

\details{
  \describe{
    \item{\code{pntR1()}:}{a pure \R version of the (C level)
      code of \R's own \code{\link{pt}()}, additionally giving more
      flexibility (via arguments \code{use.pnorm}, \code{itrmax}, \code{errmax}
      whose defaults here have been hard-coded in \R's C code).

      This implements an improved version of the AS 243 algorithm from
      Lenth(1989); \describe{
	\item{\R's help on non-central \code{\link{pt}()} says:}{
	  \emph{This computes the lower tail only, so the upper tail suffers from
	    cancellation and a warning will be given when this is likely to be
	    significant.}}
	\item{and (in \sQuote{Note:})}{\emph{The code for non-zero
	    \code{ncp} is principally intended to be used for moderate
	    values of \code{ncp}: it will not be highly accurate,
	    especially in the tails, for large values.}}
	}
    }
    \item{\code{pntR()}:}{the \code{\link{Vectorize}()}d version of \code{pntR1()}.}

    \item{\code{pntP94()}, \code{pntP94.1()}:}{New versions of
      \code{pntR1()}, \code{pntR()}; using the  Posten (1994) algorithm.
      \code{pntP94()} is the \code{\link{Vectorize}()}d version of
      \code{pntP94.1()}.}

    \item{\code{pnt3150()}, \code{pnt3150.1()}:}{
      Simple inefficient but hopefully correct version of pntP94..()
      This is really a direct implementation of formula
      (31.50), p.532 of Johnson, Kotz and Balakrishnan (1995)
    }

    \item{\code{pntLrg()}:}{provides the \code{\link{pnorm}()}
      approximation (to the non-central \eqn{t}) from
      Abramowitz and Stegun (26.7.10), p.949; which should be employed only for
      \emph{large} \code{df} and/or \code{ncp}.}

    \item{\code{pntJW39.0()}:}{use the Jennett & Welch (1939) approximation
      see Johnson et al. (1995), p. 520, after (31.26a).  This is still
      \emph{fast} for huge \code{ncp} but has \emph{wrong} asymptotic tail
      for \eqn{|t| \to \infty}{|t| -> ∞}.  Crucially needs \eqn{b=}\code{\link{b_chi}(df)}.}
    \item{\code{pntJW39()}:}{is an improved version of \code{pntJW39.0()},
      using \eqn{1-b =}\code{b_chi(df, one.minus=TRUE)} to avoid
      cancellation when computing \eqn{1 - b^2}.}

    %% \item{\code{pntChShP94()}:}{ .. }
    %% \item{\code{pntChShP94.1()}:}{ .. }
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

  Lenth, R. V. (1989). \emph{Algorithm AS 243} ---
  Cumulative distribution function of the non-central \eqn{t} distribution,
  \emph{JRSS C (Applied Statistics)} \bold{38}, 185--189.

  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \cr Formula (26.7.10), p.949

  %%__NOT_YET__
  %% Chattamvelli, R. and Shanmugam, R. (1994)
  %% An enhanced algorithm for noncentral t-distribution,
  %% \emph{Journal of Statistical Computation and Simulation} \bold{49}, 77--83.
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
