%%------------ = Part 1 of  "pnchisq Approximations"  [part 2 : ==> ./pnchisqWienergerm.Rd]
\name{pnchisqAppr}
\title{(Approximate) Probabilities of Non-Central Chi-squared Distribution}
%%     ----------------------------------------------------------------
\alias{pnchisq}
\alias{pnchisqV}% with "n" (--> ./pnchisqWienergerm.Rd has {pchisqV}  w/o "n" )
\alias{pnchisqRC}
\alias{pnchisq_ss}
\alias{pnchisqAbdelAty}
\alias{pnchisqBolKuz}
\alias{pnchisqIT}
\alias{pnchisqPatnaik}
\alias{pnchisqPearson}
\alias{pnchisqSankaran_d}
\alias{pnchisqTerms}
\alias{pnchisqT93}
\alias{pnchisqT93.a}
\alias{pnchisqT93.b}
%
\alias{ss}
\alias{ss2}
\alias{ss2.}
\description{
  Compute (approximate) probabilities for the non-central chi-squared
  distribution.

  The non-central chi-squared distribution with \code{df}\eqn{= n}
  degrees of freedom and non-centrality parameter \code{ncp}
  \eqn{= \lambda} has density
  \deqn{
    f(x) = f_{n,\lambda}(x) = e^{-\lambda / 2}
      \sum_{r=0}^\infty \frac{(\lambda/2)^r}{r!}\, f_{n + 2r}(x)}{f(x) = exp(-\lambda/2) SUM_{r=0}^\infty ((\lambda/2)^r / r!) dchisq(x, df + 2r)
  }
  for \eqn{x \ge 0}; for more, see \R's help page for \code{\link{pchisq}}.
  % ~/R/D/r-devel/R/src/library/stats/man/Chisquare.Rd

  \itemize{
    \item \R's own historical and current versions, but with more tuning parameters;
  }
  Historical relatively simple approximations listed in Johnson, Kotz, and Balakrishnan (1995):
  \itemize{
    \item Patnaik(1949)'s approximation to the non-central via central
    chi-squared.  Is also the formula \eqn{26.4.27} in Abramowitz & Stegun, p.942.
    Johnson et al mention that the approximation error is \emph{
      \eqn{O(1/\sqrt(\lambda))}{O(1/sqrt(ncp))} for \eqn{\lambda \to \infty}{ncp -> Inf}}.

    \item Pearson(1959) is using 3 moments instead of 2 as Patnaik (to
    approximate via a central chi-squared), and therefore better than
    Patnaik for the right tail; further (in Johnson et al.), the
    approximation error is \emph{\eqn{O(1/\lambda)}{O(1/ncp)} for \eqn{\lambda \to \infty}{ncp -> Inf}}.

    \item Abdel-Aty(1954)'s \dQuote{first approximation} based on
    Wilson-Hilferty via Gaussian (\code{\link{pnorm}}) probabilities, is
    partly \emph{wrongly} cited in Johnson et al., p.463, eq.\eqn{(29.61a)}.

    \item Bol'shev and Kuznetzov (1963) concentrate on the case of
    \bold{small} \code{ncp} \eqn{\lambda} and provide an \dQuote{approximation} via
    \emph{central} chi-squared with the same degrees of freedom \code{df},
    but a modified \code{q} (\sQuote{x}); the approximation has error
    \eqn{O(\lambda^3)} for \eqn{\lambda \to 0}{\lambda -> 0} and is from
    Johnson et al., p.465, eq.\eqn{(29.62)} and \eqn{(29.63)}.

    \item Sankaran(1959, 1963) proposes several further approximations base
    on Gaussian probabilities, according to Johnson
    et al., p.463. \code{pnchisqSankaran_d()} implements its formula \eqn{(29.61d)}.
  }

  \describe{
    \item{\code{pnchisq()}:}{an R implementation of \R's own C \code{pnchisq_raw()},
      but almost only up to Feb.27, 2004, long before the \code{log.p=TRUE}
      addition there, including \emph{logspace arithmetic} in April 2014,
      its finish on 2015-09-01.  Currently for historical reference only.}
    %% but we could add the log space computations from C to the R code
    \item{\code{pnchisqV()}:}{a \code{\link{Vectorize}()}d \code{\link{pnchisq}}.}

    \item{\code{pnchisqRC()}:}{\R's C implementation as of Aug.2019; but
      with many more options.
      Currently extreme cases tend to hang on Winbuilder (?) % FIXME?
    }

    \item{\code{pnchisqIT}:}{ .... }
    \item{\code{pnchisqTerms}:}{ .... }
    %
    \item{\code{pnchisqT93}:}{pure \R implementations of approximations when
      both \code{q} and \code{ncp} are large, by Temme(1993), from Johnson
      et al., p.467, formulas \eqn{(29.71 a)}, and \eqn{(29.71 b)}, using
      auxiliary functions \code{pnchisqT93a()} and \code{pnchisqT93b()}
      respectively, with adapted formulas for the \code{log.p=TRUE} cases.}

    \item{\code{pnchisq_ss()}:}{ .... }
    %
    \item{\code{ss}:}{ .... }
    \item{\code{ss2}:}{ .... }
    \item{\code{ss2.}:}{ .... }
  }% end{descr}
}
\usage{
%% R code in  ../R/pnchisq.R :
%%                ----------
pnchisq          (q, df, ncp = 0, lower.tail = TRUE, %% TODO?!? log.p = FALSE,
                  cutOffncp = 80, itSimple = 110, errmax = 1e-12, reltol = 1e-11,
                  maxit = 10* 10000, verbose = 0, xLrg.sigma = 5)
pnchisqV(x, \dots, verbose = 0)

pnchisqRC        (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE,
                  no2nd.call = FALSE,
                  cutOffncp = 80, small.ncp.logspace = small.ncp.logspaceR2015,
                  itSimple = 110, errmax = 1e-12,
                  reltol = 8 * .Machine$double.eps, epsS = reltol/2, maxit = 1e6,
                  verbose = FALSE)
pnchisqAbdelAty  (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
pnchisqBolKuz    (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
pnchisqPatnaik   (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
pnchisqPearson   (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
pnchisqSankaran_d(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
pnchisq_ss       (x, df, ncp = 0, lower.tail = TRUE, log.p = FALSE, i.max = 10000)
pnchisqTerms     (x, df, ncp,     lower.tail = TRUE, i.max = 1000)

pnchisqT93  (q, df, ncp, lower.tail = TRUE, log.p = FALSE, use.a = q > ncp)
pnchisqT93.a(q, df, ncp, lower.tail = TRUE, log.p = FALSE)
pnchisqT93.b(q, df, ncp, lower.tail = TRUE, log.p = FALSE)

ss   (x, df, ncp, i.max = 10000, useLv = !(expMin < -lambda && 1/lambda < expMax))
ss2  (x, df, ncp, i.max = 10000, eps = .Machine$double.eps)
ss2. (q, df, ncp = 0, errmax = 1e-12, reltol = 2 * .Machine$double.eps,
      maxit = 1e+05, eps = reltol, verbose = FALSE)
}
\arguments{%--> ~/R/D/r-devel/R/src/library/stats/man/Chisquare.Rd :
%% ===>  pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  \item{x}{numeric vector (of \sQuote{quantiles}, i.e., abscissa values).}
  \item{q}{number ( \sQuote{quantile}, i.e., abscissa value.)}
  \item{df}{degrees of freedom \eqn{> 0}, maybe non-integer.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta}; ....}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{pchisq}()}.}
%%% TODO: \item{method}{a string specifying the approximation method to be used.}
%%% -----  INSTEAD of all these functions ??
  \item{i.max}{number of terms in evaluation ...}
  \item{use.a}{\code{\link{logical}} vector for Temme \code{pnchisqT93*()}
    formulas, indicating to use formula \sQuote{a} over \sQuote{b}.  The
    default is as recommended in the references, but they did not take into
    account \code{log.p = TRUE} situations.}
  \item{cutOffncp}{a positive number, the cutoff value for \code{ncp}...}%% FIXME
  \item{itSimple}{ ...}
  \item{errmax}{absolute error tolerance.}% ...
  \item{reltol}{convergence tolerance for \emph{rel}ative error.}%... TODO
  \item{maxit}{maximal number of iterations.}% ...
  \item{xLrg.sigma}{positive number ...}% ... TODO
  \item{no2nd.call}{logical indicating if a 2nd call is made to the
    internal function ....}%TODO
  \item{small.ncp.logspace}{logical vector or \code{\link{function}},
    indicating if the logspace computations for \dQuote{small} \code{ncp}
    (defined to fulfill \code{ncp < cutOffncp} !).}
  \item{epsS}{small positive number, the convergence tolerance of the
    \sQuote{simple} iterations...}% ... TODO
  \item{verbose}{logical or integer specifying if or how much the algorithm
    progress should be monitored.}
  \item{\dots}{further arguments passed from \code{pnchisqV()} to \code{pnchisq()}.}
  \item{useLv}{\code{\link{logical}} indicating if logarithmic scale should
    be used for \eqn{\lambda} computations.}
  \item{eps}{convergence tolerance, a positive number.}
}
\details{
  \describe{
    \item{\code{pnchisq_ss()}}{uses \code{si <- ss(x, df, ..)} to get the series terms,
      and returns \code{2*\link{dchisq}(x, df = df +2) * sum(si$s)}.}

    \item{\code{ss()}}{computes the terms needed for the expansion used in
      \code{pnchisq_ss()}.}

    \item{\code{ss2()}}{computes some simple \dQuote{statistics} about \code{ss(..)}.}
    %% do we need it?  currently it computes ss() and trashes most of it ..
  }
}
\value{
  \describe{
    \item{\code{ss()}}{returns a list with 3 components}
  }
  \item{s}{the series}
  \item{i1}{location (in \code{s[]}) of the first change from 0 to positive.}
  \item{max}{(first) location of the maximal value in the series (i.e.,
    \code{\link{which.max}(s)}).}

  %% ss2 :  \item{i2}{the length of \code{s[]}, \code{\link{length}(s)}.}

}
\references{
  Johnson, N.L., Kotz, S. and Balakrishnan, N. (1995)
  Continuous Univariate Distributions Vol~2, 2nd ed.; Wiley;
  % \doi{10...............}
  \cr Chapter 29 \emph{Noncentral \eqn{\chi^2}-Distributions};
  notably Section \emph{8  Approximations}, p.461 ff.

  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun}
}
\author{Martin Maechler, from May 1999;  starting from a post to the S-news
  mailing list by Ranjan Maitra (@ math.umbc.edu) who showed a version of
  our \code{pchisqAppr.0()} thanking Jim Stapleton for providing it.
}
\seealso{
 \code{\link{pchisq}} and the wienergerm approximations for it:
 \code{\link{pchisqW}()} etc.

 \code{\link{r_pois}()} and its plot function, for an aspect of the series
 approximations we use in \code{pnchisq_ss()}.
}
\examples{
## set of quantiles to use :
qq <- c(.001, .005, .01, .05, (1:9)/10, 2^seq(0, 10, by= 0.5))
## Take "all interesting" pchisq-approximation from our pkg :
pkg <- "package:DPQ"
pnchNms <- c(paste0("pchisq", c("V", "W", "W.", "W.R")),
             ls(pkg, pattern = "^pnchisq"))
pnchNms <- pnchNms[!grepl("Terms$", pnchNms)]
pnchF <- sapply(pnchNms, get, envir = as.environment(pkg))
str(pnchF)
ncps <- c(0, 1/8, 1/2)
pnchR <- as.list(setNames(ncps, paste("ncp",ncps, sep="=")))
for(i.n in seq_along(ncps)) {
  ncp <- ncps[i.n]
  pnF <- if(ncp == 0) pnchF[!grepl("chisqT93", pnchNms)] else pnchF
  pnchR[[i.n]] <- sapply(pnF, function(F)
            Vectorize(F, names(formals(F))[[1]])(qq, df = 3, ncp=ncp))
}
str(pnchR, max=2)%% FIXME: *look* more; do plots , .... and/or move this
		 %% ------  to ../tests/chisq-nonc-ex.R

## A case where the non-central P[] should be improved :
## First, the central P[] which is close to exact -- choosing df=2 allows
## truly exact values: chi^2 = Exp(1) !
opal <- palette()
palette(c("black", "red", "green3", "blue", "cyan", "magenta", "gold3", "gray44"))
cR  <- curve(pchisq   (x, df=2,        lower.tail=FALSE, log.p=TRUE), 0, 4000, n=2001)
cRC <- curve(pnchisqRC(x, df=2, ncp=0, lower.tail=FALSE, log.p=TRUE),
             add=TRUE, col=adjustcolor(2,1/2), lwd=3, lty=2, n=2001)
cR0 <- curve(pchisq   (x, df=2, ncp=0, lower.tail=FALSE, log.p=TRUE),
             add=TRUE, col=adjustcolor(3,1/2), lwd=4,        n=2001)
## smart "named list" constructur :
list_ <- function(...)
   `names<-`(list(...), vapply(sys.call()[-1L], as.character, ""))
JKBfn <-list_(pnchisqPatnaik,
              pnchisqPearson,
              pnchisqAbdelAty,
              pnchisqBolKuz,
              pnchisqSankaran_d)
cl. <- setNames(adjustcolor(3+seq_along(JKBfn), 1/2), names(JKBfn))
lw. <- setNames(2+seq_along(JKBfn),                   names(JKBfn))
cR.JKB <- sapply(names(JKBfn), function(nmf) {
  curve(JKBfn[[nmf]](x, df=2, ncp=0, lower.tail=FALSE, log.p=TRUE),
        add=TRUE, col=cl.[[nmf]], lwd=lw.[[nmf]], lty=lw.[[nmf]], n=2001)
})
legend("bottomleft", c("pchisq", "pchisq.ncp=0", "pnchisqRC", names(JKBfn)),
       col=c(palette()[1], adjustcolor(2:3,1/2), cl.),
       lwd=c(1,3,4, lw.), lty=c(1,2,1, lw.))
palette(opal)# revert

all.equal(cRC, cR0, tol = 1e-15) # TRUE [for now]
## the problematic "jump" :
as.data.frame(cRC)[744:750,]
if(.Platform$OS.type == "unix")
  ## verbose=TRUE  may reveal which branches of the algorithm are taken:
  pnchisqRC(1500, df=2, ncp=0, lower.tail=FALSE, log.p=TRUE, verbose=TRUE) #
  ## |-->  -Inf currently

## The *two*  principal cases (both lower.tail = {TRUE,FALSE} !), where
##  "2nd call"  happens *and* is currently beneficial :
dfs <- c(1:2, 5, 10, 20)
pL. <- pnchisqRC(.00001, df=dfs, ncp=0, log.p=TRUE, lower.tail=FALSE, verbose = TRUE)
pR. <- pnchisqRC(   100, df=dfs, ncp=0, log.p=TRUE,                   verbose = TRUE)
## R's own non-central version (specifying 'ncp'):
pL0 <- pchisq   (.00001, df=dfs, ncp=0, log.p=TRUE, lower.tail=FALSE)
pR0 <- pchisq   (   100, df=dfs, ncp=0, log.p=TRUE)
## R's *central* version, i.e., *not* specifying 'ncp' :
pL  <- pchisq   (.00001, df=dfs,        log.p=TRUE, lower.tail=FALSE)
pR  <- pchisq   (   100, df=dfs,        log.p=TRUE)
cbind(pL., pL, relEc = signif(1-pL./pL, 3), relE0 = signif(1-pL./pL0, 3))
cbind(pR., pR, relEc = signif(1-pR./pR, 3), relE0 = signif(1-pR./pR0, 3))
}
\keyword{distribution}
