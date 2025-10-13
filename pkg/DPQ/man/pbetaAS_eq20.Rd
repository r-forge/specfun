\name{pbetaAS_eq20}
\title{\code{pbeta()} Approximations of Abramowitz & Stegun, \eqn{26.5.\{20,21\}}}
\alias{pbetaAS_eq20}
\alias{pbetaAS_eq21}
\alias{pbetaNorm2}
% "hidden" lowlevel
\alias{.pbeta.eq20}
\alias{.pbetaSeq20}
\alias{.pbeta.eq21}
\alias{.pbetaSeq21}
\description{
  Compute (simple) \code{pbeta()} approximations; Abramowitz & Stegun,
  eq. 26.5.20 and 26.5.21.
  The first, eq. 26.5.20, approximates via a \eqn{\chi^2} distribution,
  i.e., calls \code{\link{pchisq}()} with appropriate arguments,
  where the 2nd, eq. 26.5.21., approximates via the normal (aka
  \dQuote{Gaussian}) distribution, i.e., calls \code{\link{pnorm}()}.

  In addition, rather for \dQuote{didactical} reasons or simple
  comparisons, we also provide \code{pbetaNorm2()} which simply uses
  the normal approximation with matching moments, i.e., uses
  \eqn{Phi^{-1}((x - \mu)/\sigma)} where \eqn{\mu = \mu(\alpha,\beta)} is
  the mean, and   \eqn{\sigma = \sigma(\alpha,\beta)} is the standard
  deviation of the \eqn{B(\alpha,\beta)} Beta distribution, where
  \eqn{\alpha=}\code{shape1}, \eqn{\beta=}\code{shape2},
  \deqn{\mu(\alpha,\beta) = \frac{\alpha}{\alpha+\beta},\ \ \mathrm{and}}{%
        mu(a, b) = a/(a+b),  and}
  \deqn{\sigma(\alpha,\beta) = \frac{\sqrt{\alpha\beta},(\alpha+\beta)\sqrt{\alpha+\beta+1}}.}{%
         sigma(a,b) = sqrt(ab)/(a+b)/sqrt(a+b+1).}
  If the \dQuote{validity} conditions are not fulfilled, the functions will
  try to \dQuote{swap tails}, i.e., using the (anti-) symmetries of the
  Beta distribution, viz. computing
  \eqn{1 - F(1-x, b, a)} instead of \eqn{F(x, a, b)}.
}
\usage{
pbetaAS_eq20(q, shape1, shape2, lower.tail = TRUE, log.p = FALSE, warnBad = TRUE)
pbetaAS_eq21(q, shape1, shape2, lower.tail = TRUE, log.p = FALSE, warnBad = TRUE)
pbetaNorm2  (q, shape1, shape2, lower.tail = TRUE, log.p = FALSE)

## lower-level -- direct formula without range checks :
.pbeta.eq20(q, a,b, lower.tail = TRUE, log.p = FALSE, ._1_q = .5 - q + .5)
.pbetaSeq20(q, a,b, lower.tail = TRUE, log.p = FALSE)
.pbeta.eq21(q, a,b, lower.tail = TRUE, log.p = FALSE, ._1_q = .5 - q + .5)
.pbetaSeq21(q, a,b, lower.tail = TRUE, log.p = FALSE)

}
\arguments{
  \item{q, shape1, shape2}{non-negative vectors (recycled to common
    length), \code{q} in \eqn{[0,1]}, see \code{\link{pbeta}}.}
  \item{a, b}{the same as \code{shape1} and \code{shape2}, respectively;
    for convenience in (\R) formulas only.}
  \item{lower.tail}{indicating if \eqn{F(q; *)} should be returned or
    the upper tail probability \eqn{1 - F(q)}.}
  \item{log.p}{logical indicating if \code{\link{log}(<pr>)},
    log-probabilities should be returned instead of (default) the
    probabilities \code{<pr>}.}
  \item{warnBad}{logical indicating if a \code{\link{warning}} should be
    signalled when for some \code{(q, shape1, shape2, lower.tail)}
    combinations, A.&S. indicate the approximation to be bad, i.e.,
    potentially quite inaccurate.}
  \item{._1_q}{\eqn{ = 1 - q}, computed numerically carefully (potentially
    with less cancellation error), or pre-computed.}
}
\value{
  a numeric vector, the length of e.g., \code{length(q + shape1 + shape2)},
  i.e., recycled to common length.
  With the corresponding \code{pbeta()} approximations, i.e., probabilities
  (in \eqn{[0,1]}) or log-probabilities in \eqn{[-\infty, 0]}.
}
\references{
  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain;
  formulas (26.5.20) & (26.5.21), p.945.

  Beta Distribution (in Wikipedia)
  \url{https://en.wikipedia.org/wiki/Beta_distribution}
}
\author{Martin Maechler}
%% \note{
%% }
\seealso{
  \code{\link{pbeta}}; our \code{\link{pbetaRv1}}, and approximations to
  the non-central beta \code{\link{pnbetaAppr2}}.
}
\examples{
str(alf <- 2^seq(1,15, by=1/8)) # 2 1.18 ... 32768
bet <- seq_along(alf) # 1 2 .. 113
% ------- TODO: use two smallish computation and plot() functions ---%
E <- alf/(alf + bet) # = E[X]
pb   <- pbeta       (E, alf, bet)
pb20 <- pbetaAS_eq20(E, alf, bet) # warning .. 112 potentially bad
pb21 <- pbetaAS_eq21(E, alf, bet)
pbN2 <- pbetaNorm2  (E, alf, bet); summary(pbN2)# all == 0.5 as q = E[ Beta(a,b) ]
pbM <- cbind(pb20, pb21, pbN2)
tit <- quote(list(pbeta[...](E, alpha, beta) ~~" ", alpha == 2^{1~".."~15},
                  beta==1:113, E == frac(alpha,alpha+beta)))
matplot(alf, cbind(pb, pbM), type = "l", lwd = 1:4, log = "x",
        xlab = quote(alpha ~~ "[log]"), xaxt = "n", main = tit)
sfsmisc::eaxis(1, sub10=2)
legend("topright", ncol=2, bty="n", lty = 1:4, col = adjustcolor(1:4, 3/4), lwd = 1:4,
       expression(pbeta(), pbetaAS_eq20(), pbetaAS_eq21(), pbetaNorm2()))
stopifnot(apply(pbM, 2, all.equal, target = pb, tolerance = 0.03))
## relative errors, compared to R's pbeta()
round(cbind(sort(apply(pbM, 2, sfsmisc::relErr, target = pb))), 5)
## pb21 0.00013 << good here
## pb20 0.00356 << "ok" (got warning about "bad")
## pbN2 0.02742 << (constant)

## 2.  below E[]: mu - 2 s =======================================================
## -----------------------  --> use  xmin :
sdB <- sqrt(alf * bet)/(alf + bet)/sqrt(alf + bet + 1)
xmin <- unique(pmax(0, E - 2*sdB))
xmax <- unique(pmin(1, E + 2*sdB))
##
pb   <- pbeta       (xmin, alf, bet)
pb20 <- pbetaAS_eq20(xmin, alf, bet) # warning .. 104 potentially bad
## --> analysis shows the 'swap' logical to be clearly sub optimal.  ==> functions .pbeta.eq20() etc
pb21 <- pbetaAS_eq21(xmin, alf, bet)
pbN2 <- pbetaNorm2  (xmin, alf, bet); summary(pbN2)
pbM <- cbind(pb20, pb21, pbN2)
round(cbind(sort(apply(pbM, 2, sfsmisc::relErr, target = pb))), 5)
## pb21 0.00594
## pbN2 0.27069
## pb20 0.68668

tit <- quote(list(pbeta[...](xmin, alpha, beta) ~~" ", alpha == 2^{1~".."~15},
                  beta==1:113, xmin == (mu - 2*sigma)(alpha,beta)))
matplot(alf, cbind(pb, pbM), type = "l", lwd = 1:4, log = "x",
        xlab = quote(alpha ~~ "[log]"), xaxt = "n", main = tit)
sfsmisc::eaxis(1, sub10=2)
legend("topright", ncol=2, bty="n", lty = 1:4, lwd = 1:4, col = 1:4, # col = adjustcolor(1:4, 3/4),
       legend = expression(pbeta(), pbetaAS_eq20(), pbetaAS_eq21(), pbetaNorm2()))
##
pb20. <- .pbeta.eq20(xmin, alf,bet)
pb20S <- .pbetaSeq20(xmin, alf,bet)
pb21. <- .pbeta.eq21(xmin, alf,bet)
pb21S <- .pbetaSeq21(xmin, alf,bet)
matlines(alf, cbind(pb20., pb20S, pb21., pb21S), lwd = 4, lty = 1, col = adjustcolor(2:5, 1/4))
legend("right", ncol=2, bty="n", lty = 1, lwd = 4, col = adjustcolor(2:5, 1/4),
       legend = paste0(".pbeta", c(".eq20", "Seq20", ".eq21", "Seq21"), "()"))

## 3.  above E[]: mu + 2 s =======================================================
## -----------------------  --> use  xmax
##
pb   <- pbeta       (xmax, alf, bet)
pb20 <- pbetaAS_eq20(xmax, alf, bet) # warning .. 104 potentially bad
## --> analysis shows the 'swap' logical to be clearly sub optimal.  ==> functions .pbeta.eq20() etc
pb21 <- pbetaAS_eq21(xmax, alf, bet)
pbN2 <- pbetaNorm2  (xmax, alf, bet); summary(pbN2)
pbM <- cbind(pb20, pb21, pbN2)
round(cbind(sort(apply(pbM, 2, sfsmisc::relErr, target = pb))), 5)
## pb21 0.00009
## pb20 0.00559
## pbN2 0.00610
tit <- quote(list(pbeta[...](xmax, alpha, beta) ~~" ", alpha == 2^{1~".."~15},
                  beta==1:113, xmax == (mu + 2*sigma)(alpha,beta)))
matplot(alf, cbind(pb, pbM), type = "l", lwd = 1:4, log = "x",
        xlab = quote(alpha ~~ "[log]"), xaxt = "n", main = tit)
sfsmisc::eaxis(1, sub10=2)
legend("topright", ncol=2, bty="n", lty = 1:4, lwd = 1:4, col = 1:4, # col = adjustcolor(1:4, 3/4),
       legend = expression(pbeta(), pbetaAS_eq20(), pbetaAS_eq21(), pbetaNorm2()))
pb20. <- .pbeta.eq20(xmax, alf,bet)
pb20S <- .pbetaSeq20(xmax, alf,bet)
pb21. <- .pbeta.eq21(xmax, alf,bet)
pb21S <- .pbetaSeq21(xmax, alf,bet)
matlines(alf, cbind(pb20., pb20S, pb21., pb21S), lwd = 4, lty = 1, col = adjustcolor(2:5, 1/4))
legend("right", ncol=2, bty="n", lty = 1, lwd = 4, col = adjustcolor(2:5, 1/4),
       legend = paste0(".pbeta", c(".eq20", "Seq20", ".eq21", "Seq21"), "()"))

}
\keyword{math}
\keyword{distribution}
