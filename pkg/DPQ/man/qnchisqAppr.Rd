\name{qnchisqAppr}
\alias{qnchisqAppr}% < even though have no such fn.name
\title{Compute Approximate Quantiles of Noncentral Chi-Squared Distribution}
%%     --------------------------------------------------------------------
\alias{qchisqAppr.0}
\alias{qchisqAppr.1}
\alias{qchisqAppr.2}
\alias{qchisqAppr.3}
\alias{qchisqApprCF1}
\alias{qchisqApprCF2}
%-
\alias{qchisqCappr.2}
\alias{qchisqN}
%-
\alias{qnchisqAbdelAty}
\alias{qnchisqPatnaik}
\alias{qnchisqPearson}
\alias{qnchisqSankaran_d}
\description{%% R code:    ../R/qnchisq.R
  %            %           --------------
  Compute quantiles (inverse distribution values) for the \emph{non-central}
  chi-squared distribution.

  ....... using Johnson,Kotz,  and other approximations ..............
}
\usage{
qchisqAppr.0 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqAppr.1 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqAppr.2 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqAppr.3 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqApprCF1(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqApprCF2(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)

qchisqCappr.2 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqN       (p, df, ncp = 0, qIni = qchisqAppr.0, \dots)

qnchisqPatnaik   (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qnchisqAbdelAty  (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qnchisqPearson   (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qnchisqSankaran_d(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
}
\arguments{%--> ~/R/D/r-devel/R/src/library/stats/man/Chisquare.Rd :
%% ===>  qchisq(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  \item{p}{vector of probabilities.}
  \item{df}{degrees of freedom \eqn{> 0}, maybe non-integer.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta}; ....}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qchisq}()}.}
%%% TODO: \item{method}{a string specifying the approximation method to be used.}
%%% -----  INSTEAD of all these functions ??
  \item{qIni}{a \code{\link{function}} that computes an approximate
    noncentral chi-squared quantile as starting value \code{x0} for the
    Newton algorithm \code{\link{newton}()}.}
  \item{\dots}{further arguments to \code{\link{newton}()}, notably
    \code{eps} or \code{maxiter}.}
}
\details{
  \describe{
    \item{qchisqAppr.0():}{ ...TODO... }
    \item{qchisqAppr.1():}{ ...TODO... }
    \item{qchisqAppr.2():}{ ...TODO... }
    \item{qchisqAppr.3():}{ ...TODO... }
    \item{qchisqApprCF1():}{ ...TODO... }
    \item{qchisqApprCF2():}{ ...TODO... }

    \item{qchisqCappr.2():}{ ...TODO... }
    \item{qchisqN():}{Uses Newton iterations with \code{\link{pchisq}()}
      and \code{\link{dchisq}()} to determine \code{\link{qchisq}(.)} values.}

    \item{qnchisqPatnaik():}{ ...TODO... }
    \item{qnchisqAbdelAty():}{ ...TODO... }
    \item{qnchisqPearson():}{ ...TODO... }
    \item{qnchisqSankaran_d():}{ ...TODO... }

  }
}
\value{
  ...
}
%% \references{
%%
%% }
\author{Martin Maechler, from May 1999;  starting from a post to the S-news
  mailing list by Ranjan Maitra (@ math.umbc.edu) who showed a version of
  our \code{qchisqAppr.0()} thanking Jim Stapleton for providing it.
}
\seealso{
 \code{\link{qchisq}}.
}
\examples{
 pp <- c(.001, .005, .01, .05, (1:9)/10, .95, .99, .995, .999)
 pkg <- "package:DPQ"
 qnchNms <- c(paste0("qchisqAppr.",0:3), paste0("qchisqApprCF",1:2),
              "qchisqN", "qchisqCappr.2", ls(pkg, pattern = "^qnchisq"))
 qnchF <- sapply(qnchNms, get, envir = as.environment(pkg))
 for(ncp in c(0, 1/8, 1/2)) {
   cat("\n~~~~~~~~~~~~~\nncp: ", ncp,"\n=======\n")
   print(sapply(qnchF, function(F) Vectorize(F, "p")(pp, df = 3, ncp=ncp)))
 }
}
\keyword{distribution}
