\name{pnt}
\title{Non-central t Probability Distribution - Algorithms and Approximations}
\alias{pntR}
\alias{pntR1}
\alias{pntP94}
\alias{pntP94.1}
\alias{pntLrg}
\alias{pntJW39}
\alias{pntJW39.0}
\alias{pnt3150}
\alias{pnt3150.1}
\alias{pntChShP94}
\alias{pntChShP94.1}
\alias{pntGST23_T6}
\alias{pntGST23_T6.1}
\alias{pntGST23_1}
\alias{pntVW13}
\description{
 Compute different approximations for the non-central t-Distribution
 cumulative probability distribution function.
}
\usage{% --> ~/R/D/r-devel/R/src/library/stats/man/TDist.Rd
%% source in ../R/t-nonc-fn.R
pntR      (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
           use.pnorm = (df > 4e5 ||
                        ncp^2 > 2*log(2)*1021), # .Machine$double.min.exp = -1022
                                          itrmax = 1000, errmax = 1e-12, verbose = TRUE)
pntR1     (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
           use.pnorm = (df > 4e5 ||
                        ncp^2 > 2*log(2)*1021),
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

pntVW13 (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
           keepS = FALSE, verbose = FALSE)

pntGST23_T6  (t, df, ncp, lower.tail = TRUE, log.p = FALSE,
              y1.tol = 1e-8, Mterms = 20, alt = FALSE, verbose = TRUE)
pntGST23_T6.1(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
              y1.tol = 1e-8, Mterms = 20, alt = FALSE, verbose = TRUE)

## *Non*-asymptotic, (at least partly much) better version of R's Lenth(1998) algorithm
pntGST23_1(t, df, ncp, lower.tail = TRUE, log.p = FALSE,
           j0max = 1e4, # for now
           IxpqFUN = Ixpq,
           alt = FALSE, verbose = TRUE, \dots)
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
    is suboptimal).}

  \item{itrmax}{number of iterations / terms.}% .. say more ?
  \item{errmax}{convergence bound for the iterations.}%% could say more

  \item{verbose}{\code{\link{logical}} or integer determining the amount of
    diagnostic print out to the console.}
  \item{M}{positive integer specifying the number of terms to use in the series.}

  \item{keepS}{\code{\link{logical}} indicating if the function should
    return a \code{\link{list}} with component \code{cdf} and other
    informational elements, or just the CDF values directly (by default).}

  \item{y1.tol}{positive tolerance for warning if \eqn{y:= t^2/(t^2 + df)}
    is too close to 1 (as the formulas use \eqn{1/(1-y)}).}
  \item{Mterms}{number of summation terms for \code{pntGST23_T6()}.}%% could say more

  \item{j0max}{\emph{experimental}: large integer limiting the summation
    terms in \code{pntGST23_1()} .}

  \item{IxpqFUN}{the (scaled) incomplete beta function \eqn{I_x(p,q)} to be
    used;  currently, it defaults to the \code{Ixpq} function derived from
    Nico Temme's Maple code for \dQuote{Table 1} in Gil et al. (2023).}

  \item{alt}{\code{\link{logical}} specifying if and how log-scale should
    be used.  \bold{Experimental} and not-yet-tested.}

  \item{\dots}{further arguments passed to \code{IxpqFUN()}.}
}

\details{
  \describe{
    \item{\code{pntR1()}:}{a pure \R version of the (C level)
      code of \R's own \code{\link{pt}()}, additionally giving more
      flexibility (via arguments \code{use.pnorm}, \code{itrmax}, \code{errmax}
      whose defaults here have been hard-coded in \R's C code called by \code{\link{pt}()}).

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
      \code{pntR1()}, \code{pntR()}; using the Posten (1994) algorithm.
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
    \item{\code{pntGST23_T6()}:}{(and \code{pntGST23_T6.1()} for
      informational purposes only) use the Gil et al.(2023)'s
      approximation of their Theorem 6.}

    \item{\code{pntGST23_1()}:}{implements Gil et al.(2023)'s direct
      \code{\link{pbeta}()} based formula (1), which is very close to
      Lenth's algorithm.}

    \item{\code{pntVW13()}:}{use MM's \R translation of Viktor
      \enc{Witkowský}{Witkowsky} (2013)'s matlab implementation.}
  }
}

\value{
  a number for \code{pntJKBf1()} and \code{.pntJKBch1()}.

  a numeric vector of the same length as the maximum of the lengths of
  \code{x, df, ncp} for  \code{pntJKBf()} and \code{.pntJKBch()}.
}
\references{
  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995)
  \emph{Continuous Univariate Distributions Vol~2}, 2nd ed.; Wiley;
  chapter 31, Section \emph{5  Distribution Function}, p.514 ff

  Lenth, R. V. (1989). \emph{Algorithm AS 243} ---
  Cumulative distribution function of the non-central \eqn{t} distribution,
  \emph{JRSS C (Applied Statistics)} \bold{38}, 185--189.

  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover;
  formula (26.7.10), p.949

  Posten, Harry O. (1994)
  A new algorithm for the noncentral t distribution function,
  \emph{Journal of Statistical Computation and Simulation} \bold{51}, 79--87;
  \doi{10.1080/00949659408811623}.

  -- not yet implemented -- \cr %%__NOT_YET__
  Chattamvelli, R. and Shanmugam, R. (1994)
  An enhanced algorithm for noncentral t-distribution,
  \emph{Journal of Statistical Computation and Simulation} \bold{49}, 77--83.
  \doi{10.1080/00949659408811561}

  -- not yet implemented -- \cr %%__NOT_YET__
  Akahira, Masafumi. (1995).
  A higher order approximation to a percentage point of the noncentral t distribution,
  \emph{Communications in Statistics - Simulation and Computation} \bold{24}:3, 595--605;
  \doi{10.1080/03610919508813261}

  %% Comparing Johnson et al (1995) / Chattamvelli & Shanmugam (1994) / Akahira (1995), .....
  Michael Perakis and Evdokia Xekalaki (2003)
  On a Comparison of the Efficacy of Various Approximations of the Critical Values for Tests on the Process Capability Indices CPL, CPU, and Cpk,
  \emph{Communications in Statistics - Simulation and Computation} \bold{32}, 1249--1264;
  \doi{10.1081/SAC-120023888}

  \enc{Witkovský}{Witkovsky}, Viktor (2013)
  A Note on Computing Extreme Tail Probabilities of the Noncentral T Distribution with Large
  Noncentrality Parameter,
  \emph{Acta Universitatis Palackianae Olomucensis, Facultas Rerum Naturalium, Mathematica}
  %% GST23 have  \emph{Acta Univ Palack Olomuc Fac Rerum Natur Math.}
  \bold{52}(2), 131--143. % on arXiv --> ../Misc/pnt-Gil_etal-2023/Witkovsky_2013/arXiv-1306.5294v2/

  Gil A., Segura J., and Temme N.M. (2023)
  New asymptotic representations of the noncentral t-distribution,
  \emph{Stud Appl Math.} \bold{151}, 857--882; \doi{10.1111/sapm.12609} ;
  acronym \dQuote{GST23}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{pt}}, for \R's version of non-central t probabilities.
}
\examples{
tt <- seq(0, 10, len = 21)
ncp <- seq(0, 6, len = 31)
pt3R   <- outer(tt, ncp, pt, , df = 3)
pt3JKB <- outer(tt, ncp, pntR, df = 3)# currently verbose
stopifnot(all.equal(pt3R, pt3JKB, tolerance = 4e-15))# 64-bit Lnx: 2.78e-16


## Gil et al.(2023) -- Table 1 p.869
str(GST23_tab1 <- read.table(header=TRUE, text = "
 x     pnt_x_delta              Rel.accuracy   l_y   j_max
 5     0.7890745035061528e-20    0.20e-13    0.29178   254
 8     0.1902963697413609e-07    0.40e-12    0.13863   294
11     0.4649258368179092e-03    0.12e-09    0.07845   310
14     0.2912746016055676e-01    0.11e-07    0.04993   317
17     0.1858422833307925e-00    0.41e-06    0.03441   321
20     0.4434882973203470e-00    0.82e-05    0.02510   323"))

x1 <- c(5,8,11,14,17,20)
(p1  <- pt  (x1, df=10.3, ncp=20))
(p1R <- pntR(x1, df=10.3, ncp=20)) # verbose=TRUE  is default
all.equal(p1, p1R, tolerance=0) # 4.355452e-15 {on x86_64} as have *no* LDOUBLE on R level
stopifnot(all.equal(p1, p1R))
## NB: According to Gil et al., the first value (x=5) is really wrong
## p1.23 <- .. Gil et al., Table 1:
p1.23.11 <- pntGST23_T6(x1, df=10.3, ncp=20, Mterms = 11)
p1.23.20 <- pntGST23_T6(x1, df=10.3, ncp=20, Mterms = 20, verbose=TRUE)
                                        # ==> Mterms = 11 is good only for x=5
p1.23.50 <- pntGST23_T6(x1, df=10.3, ncp=20, Mterms = 50, verbose=TRUE)

x <- 4:40 ; df <- 10.3
ncp <- 20
p1     <- pt        (x, df=df, ncp=ncp)
pG1    <- pntGST23_1(x, df=df, ncp=ncp)
pG1.bR <- pntGST23_1(x, df=df, ncp=ncp,
                     IxpqFUN = \(x, l_x=.5-x+.5, p, q) Ixpq(x,l_x, p,q))
pG1.BR <- pntGST23_1(x, df=df, ncp=ncp,
                     IxpqFUN = \(x, l_x, p, q)   pbeta(x, p,q))
cbind(x, p1, pG1, pG1.bR, pG1.BR)
all.equal(pG1, p1,     tolerance=0) # 1.034 e-12
all.equal(pG1, pG1.bR, tolerance=0) # 2.497031 e-13
all.equal(pG1, pG1.BR, tolerance=0) # 2.924698 e-13
all.equal(pG1.BR,pG1.bR,tolerance=0)# 1.68644  e-13
stopifnot(exprs = {
    all.equal(pG1, p1,     tolerance = 4e-12)
    all.equal(pG1, pG1.bR, tolerance = 1e-12)
    all.equal(pG1, pG1.BR, tolerance = 1e-12)
  })

ncp <- 40 ## is  > 37.62 = "critical" for Lenth' algorithm

### --------- pntVW13() --------------------------------------------------
## length 1 arguments:
str(rr <- pntVW13(t = 1, df = 2, ncp = 3, verbose=TRUE, keepS=TRUE))
all.equal(rr$cdf, pt(1,2,3), tol = 0)#  "Mean relative difference: 4.956769e-12"
stopifnot( all.equal(rr$cdf, pt(1,2,3)) )

str(rr <- pntVW13(t = 1:19, df = 2,    ncp = 3,    verbose=TRUE, keepS=TRUE))
str(r2 <- pntVW13(t = 1,    df = 2:20, ncp = 3,    verbose=TRUE, keepS=TRUE))
str(r3 <- pntVW13(t = 1,    df = 2:20, ncp = 3:21, verbose=TRUE, keepS=TRUE))

pt1.10.5_T <- 4.34725285650591657e-5 # Ex. 7 of Witkovsky(2013)
pt1.10.5 <- pntVW13(1, 10, 5)
all.equal(pt1.10.5_T, pt1.10.5, tol = 0)# TRUE! (Lnx Fedora 40; 2024-07-04);
			# 3.117e-16 (Macbuilder R 4.4.0, macOS Ventura 13.3.1)
stopifnot(exprs = {
    identical(rr$cdf, r1 <- pntVW13(t = 1:19, df = 2, ncp = 3))
    identical(r1[1], pntVW13(1, 2, 3))
    identical(r1[7], pntVW13(7, 2, 3))
    all.equal(pt1.10.5_T, pt1.10.5, tol = 9e-16)# NB even tol=0 (64 Lnx)
})
## However, R' pt() is only equal for the very first
cbind(t = 1:19, pntVW = r1, pt = pt(1:19, 2,3))

}
\keyword{distribution}
\keyword{math}
