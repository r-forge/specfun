%%------------ = Part 1 of  "pnchisq Approximations"  [part 2 : ==> ./pnchisqWienergerm.Rd]
\name{pnchisqAppr}
\title{Approximate Probabilities of Non-Central Chisquared Distribution}
%%     ----------------------------------------------------------------
\alias{pnchisq}
\alias{pnchisqV}% with "n" (--> ./pnchisqWienergerm.Rd has {pchisqV}  w/o "n" )
\alias{pnchisq_ss}
\alias{pnchisqAbdelAty}
\alias{pnchisqIT}
\alias{pnchisqPatnaik}
\alias{pnchisqPearson}
\alias{pnchisqSankaran_d}
\alias{pnchisqTerms}
%
\alias{ss}
\alias{ss2}
\alias{ss2.}
\description{
  Compute (approximate) probabilities for the non-central chi squared distribution.
  using Johnson, Kotz,.. ............
}
\usage{
%% R code in  ../R/pnchisq.R :
%                 ----------
pnchisq          (q, df, ncp = 0, lower.tail = TRUE, %% TODO?!? log.p = FALSE,
                  cutOffncp = 80, itSimple = 110, errmax = 1e-12, reltol = 1e-11,
                  itrmax = 10* 10000, verbose = 1, xLrg.sigma = 5)
pnchisqV(x, \dots, verbose = 0)

pnchisqPatnaik   (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
pnchisqPearson   (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
pnchisqAbdelAty  (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
pnchisqSankaran_d(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
pnchisq_ss       (x, df, ncp, i.max = 10000)
pnchisqTerms     (x, df, ncp, lower.tail = TRUE, i.max = 1000)

ss   (x, df, ncp, i.max = 10000)
ss2  (x, df, ncp, i.max = 10000, eps = .Machine$double.eps)
ss2. (q, df, ncp = 0, errmax = 1e-12, reltol = 2 * .Machine$double.eps,
      itrmax = 1e+05, eps = reltol)
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
  \item{cutOffncp}{a positive number, the cutoff value for \code{ncp}...}%% FIXME
  \item{itSimple}{ ...}
  \item{errmax}{absolute error tolerance.}% ...
  \item{reltol}{relative tolerance.}% ...
  \item{itrmax}{maximal number of iterations.}% ...
  \item{xLrg.sigma}{positive number ...}% ... TODO
  \item{verbose}{logical or integer specifying if or how much the algorithm
    progress should be monitored.}
  \item{\dots}{further arguments passed from \code{pnchisqV()} to \code{pnchisq()}.}
}
%% \details{
%% }
\value{
  ...
}
%% \references{
%%
%% }
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
 qq <- c(.001, .005, .01, .05, (1:9)/10, 2^seq(0, 10, by= 0.5))
 pkg <- "package:DPQ"
 pnchNms <- c(paste0("pchisq", c("V", "W", "W.", "W.R")),
              ls(pkg, pattern = "^pnchisq"))
 pnchNms <- pnchNms[!grepl("Terms$", pnchNms)]
 pnchF <- sapply(pnchNms, get, envir = as.environment(pkg))
 str(pnchF)
 for(ncp in c(0, 1/8, 1/2)) {
   cat("\n~~~~~~~~~~~~~\nncp: ", ncp,"\n=======\n")
   print(sapply(pnchF, function(F) Vectorize(F, names(formals(F))[[1]])(qq, df = 3, ncp=ncp)))
 }
}
\keyword{distribution}
