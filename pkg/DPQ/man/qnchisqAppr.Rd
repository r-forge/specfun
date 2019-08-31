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

## One current (2019-08) R bug: Noncentral chi-squared quantiles on *LOG SCALE*

## a)  left/lower tail : -------------------------------------------------------
qs <- 2^seq(0,11, by=1/16)
pqL <- pchisq(qs, df=5, ncp=1, log.p=TRUE)
plot(qs, -pqL, type="l", log="xy") # + expected warning on log(0) -- all fine
qpqL <- qchisq(pqL, df=5, ncp=1, log.p=TRUE) # severe overflow :
qm <- cbind(qs, pqL, qchisq=qpqL
	, qchA.0 = qchisqAppr.0 (pqL, df=5, ncp=1, log.p=TRUE)
	, qchA.1 = qchisqAppr.1 (pqL, df=5, ncp=1, log.p=TRUE)
	, qchA.2 = qchisqAppr.2 (pqL, df=5, ncp=1, log.p=TRUE)
	, qchA.3 = qchisqAppr.3 (pqL, df=5, ncp=1, log.p=TRUE)
	, qchACF1= qchisqApprCF1(pqL, df=5, ncp=1, log.p=TRUE)
	, qchACF2= qchisqApprCF2(pqL, df=5, ncp=1, log.p=TRUE)
	, qchCa.2= qchisqCappr.2(pqL, df=5, ncp=1, log.p=TRUE)
	, qnPatnaik   = qnchisqPatnaik   (pqL, df=5, ncp=1, log.p=TRUE)
	, qnAbdelAty  = qnchisqAbdelAty  (pqL, df=5, ncp=1, log.p=TRUE)
	, qnPearson   = qnchisqPearson   (pqL, df=5, ncp=1, log.p=TRUE)
	, qnSankaran_d= qnchisqSankaran_d(pqL, df=5, ncp=1, log.p=TRUE)
)

round(qm[ qs \%in\% 2^(0:11) , -2])
#=> Approximations don't overflow but are not good enough

## b)  right/upper tail , larger ncp -------------------------------------------
qS <- 2^seq(-3, 3, by=1/8)
pqLu <- pchisq(qS, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
## using "the alternative" (here is currently identical):
identical(pqLu, (pqLu.<- log1p(-pchisq(qS, df=5, ncp=100)))) # here TRUE
plot (qS, -pqLu, type="l", log="xy") # fine
qpqLu <- qchisq(pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
cbind(qS, pqLu, pqLu, qpqLu)# # severe underflow
qchMat <- cbind(qchisq = qpqLu
	, qchA.0 = qchisqAppr.0 (pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qchA.1 = qchisqAppr.1 (pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qchA.2 = qchisqAppr.2 (pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qchA.3 = qchisqAppr.3 (pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qchACF1= qchisqApprCF1(pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qchACF2= qchisqApprCF2(pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qchCa.2= qchisqCappr.2(pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qnPatnaik   = qnchisqPatnaik   (pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qnAbdelAty  = qnchisqAbdelAty  (pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qnPearson   = qnchisqPearson   (pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	, qnSankaran_d= qnchisqSankaran_d(pqLu, df=5, ncp=100, log.p=TRUE, lower.tail=FALSE)
	)
cbind(L2err <- sort(sqrt(colSums((qchMat - qS)^2))))
##--> "Sankaran_d", "CF1" and "CF2" are good here

plot (qS, qpqLu, type = "b", log="x", lwd=2)
lines(qS, qS, col="gray", lty=2, lwd=3)
top3 <- names(L2err)[1:3]
use <- c("qchisq", top3)
matlines(qS, qchMat[, use]) # 3 of the approximations are "somewhat ok"
legend("topleft", c(use,"True"), bty="n", col=c(palette()[1:4], "gray"),
                  lty = c(1:4,2), lwd = c(2,3, 1,1,1))
}
\keyword{distribution}
