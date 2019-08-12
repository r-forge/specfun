\name{qnchisqAppr}
\title{Compute Approximate Quantiles of Non-Central Chisquared Distribution}
%%     --------------------------------------------------------------------
\alias{qchisqAppr.0}
\alias{qchisqAppr.1}
\alias{qchisqAppr.2}
\alias{qchisqAppr.3}
\alias{qchisqApprCF1}
\alias{qchisqApprCF2}
%-
\alias{qchisqCappr.2}
\alias{qchisq2}
\alias{qchisqA}
%-
\alias{qnchisqAbdelAty}
\alias{qnchisqPatnaik}
\alias{qnchisqPearson}
\alias{qnchisqSankaran_d}
\description{%% R code:    ../R/qnchisq.R
  %            %           --------------
  Compute quantiles (inverse distribution values) for the non-central chi squared distribution.
  using Johnson,Kotz,.. ............
}
\usage{
qchisqAppr.0 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqAppr.1 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqAppr.2 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqAppr.3 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqApprCF1(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisqApprCF2(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)

qchisqCappr.2 (p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
qchisq2       (p, df, ncp = 0, eps = 1e-04)
qchisqA       (p, df, ncp = 0, eps = 1e-04, qIni = qchisqAppr.0)

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
  \item{eps}{....replace..with..tol...........TODO............}
  \item{qIni}{auxiliary \code{\link{function}} ...}
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
  our \code{qchisqAppr.0()} thanking Jim Stapleton for providing it.
}
\seealso{
 \code{\link{qchisq}}.
}
\examples{
 pp <- c(.001, .005, .01, .05, (1:9)/10, .95, .99, .995, .999)
 pkg <- "package:DPQ"
 qnchNms <- c(paste0("qchisqAppr.",0:3), paste0("qchisqApprCF",1:2),
              paste0("qchisq", c("Cappr.2", "2", "A")), ls(pkg, pattern = "^qnchisq"))
 qnchF <- sapply(qnchNms, get, envir = as.environment(pkg))
 for(ncp in c(0, 1/8, 1/2)) {
   cat("\n~~~~~~~~~~~~~\nncp: ", ncp,"\n=======\n")
   print(sapply(qnchF, function(F) Vectorize(F, "p")(pp, df = 3, ncp=ncp)))
 }
}
\keyword{distribution}
