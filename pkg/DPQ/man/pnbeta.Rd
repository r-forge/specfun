\name{pnbeta}
\alias{pnbetaAppr2}
\alias{pnbetaAppr2v1}
\alias{pnbetaAS310}
\title{Noncentral Beta Probabilities}
\description{
  \code{pnbetaAppr2()} and its inital version \code{pnbetaAppr2v1()}
  provide the \dQuote{approximation 2} of Chattamvelli and Shanmugam(1997)
  to the noncentral Beta probability distribution.

  \code{pnbetaAS310()} is an \R level interface to a C translation (and
  \dQuote{Rification}) of the \code{AS 310} Fortran implementation.
}
\usage{%  R code : ../R/pnbeta-approx2.R <<<<
pnbetaAppr2(x, a, b, ncp = 0, lower.tail = TRUE, log.p = FALSE)
}% R's Beta Distrib.: ~/R/D/r-devel/R/src/library/stats/man/Beta.Rd
\arguments{
  \item{x}{numeric vector (of quantiles), typically from inside \eqn{[0,1]}.}
  \item{a, b}{the shape parameters of Beta, aka as \code{shape1} and \code{shape2}.}
  \item{ncp}{non-centrality parameter.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
%% \details{

%% }
\value{
  a numeric vector of (log) probabilities of the same length as \code{x}.
}
\references{
  Chattamvelli, R., and Shanmugam, R. (1997)
  Algorithm AS 310: Computing the Non-Central Beta Distribution Function.
  \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}
  \bold{46}(1), 146--156, for \dQuote{approximation 2} notably p.154;
  \url{https://www.jstor.org/stable/2986212} or \doi{10.1111/1467-9876.00055}.
}
\author{Martin Maechler; \code{pnbetaAppr2()} in Oct 2007.}
\note{
  The authors in the reference compare AS 310 with Frick(1990) and Lenth(1987) and
  state to be better than both.  \R's current (2019) noncentral beta
  implementation builds on Lenth(1987), with some amendments though; still, \code{pnbetaAS310()}
  potentially be better, at least in certain corners of the 4-dimensional input space.
}
\seealso{
  \R's own \code{\link{pbeta}}.
}
\examples{
## Same arguments as for Table 1 (p.151) of the reference
a <- 5*rep(1:3, each=3)
aargs <- cbind(a = a, b = a,
               ncp = rep(c(54, 140, 170), 3),
               x = 1e-4*c(8640, 9000, 9560, 8686, 9000, 9000, 8787, 9000, 9220))
aargs
pnbA2 <- apply(aargs, 1, function(aa) do.call(pnbetaAppr2, as.list(aa)))
pnA310<- apply(aargs, 1, function(aa) do.call(pnbetaAS310, as.list(aa)))
aar2 <- aargs; dimnames(aar2)[[2]] <- c(paste0("shape", 1:2), "ncp", "q")
pnbR  <- apply(aar2,  1, function(aa) do.call(pbeta, as.list(aa)))
range(relD2   <- 1 - pnbA2 /pnbR)
range(relD310 <- 1 - pnA310/pnbR)
cbind(aargs, pnbA2, pnA310, pnbR,
      relD2 = signif(relD2, 3), relD310 = signif(relD310, 3)) # <------> Table 1
stopifnot(abs(relD2)   < 0.009) # max is 0.006286
stopifnot(abs(relD310) < 1e-5 ) # max is 6.3732e-6

## Arguments as for Table 2 (p.152) of the reference :
aarg2 <- cbind(a = c( 10, 10, 15, 20, 20, 20, 30, 30),
               b = c( 20, 10,  5, 10, 30, 50, 20, 40),
               ncp=c(150,120, 80,110, 65,130, 80,130),
               x = c(868,900,880,850,660,720,720,800)/1000)
pnbA2 <- apply(aarg2, 1, function(aa) do.call(pnbetaAppr2, as.list(aa)))
pnA310<- apply(aarg2, 1, function(aa) do.call(pnbetaAS310, as.list(aa)))
aar2 <- aarg2; dimnames(aar2)[[2]] <- c(paste0("shape", 1:2), "ncp", "q")
pnbR  <- apply(aar2,  1, function(aa) do.call(pbeta, as.list(aa)))
range(relD2   <- 1 - pnbA2 /pnbR)
range(relD310 <- 1 - pnA310/pnbR)
cbind(aarg2, pnbA2, pnA310, pnbR,
      relD2 = signif(relD2, 3), relD310 = signif(relD310, 3)) # <------> Table 2
stopifnot(abs(relD2  ) < 0.006) # max is 0.00412
stopifnot(abs(relD310) < 1e-5 ) # max is 5.5953e-6

## Arguments as for Table 3 (p.152) of the reference :
aarg3 <- cbind(a = c( 10, 10, 10, 15, 10, 12, 30, 35),
               b = c(  5, 10, 30, 20,  5, 17, 30, 30),
               ncp=c( 20, 54, 80,120, 55, 64,140, 20),
               x = c(644,700,780,760,795,560,800,670)/1000)
pnbA3 <- apply(aarg3, 1, function(aa) do.call(pnbetaAppr2, as.list(aa)))
pnA310<- apply(aarg3, 1, function(aa) do.call(pnbetaAS310, as.list(aa)))
aar3 <- aarg3; dimnames(aar3)[[2]] <- c(paste0("shape", 1:2), "ncp", "q")
pnbR  <- apply(aar3,  1, function(aa) do.call(pbeta, as.list(aa)))
range(relD2   <- 1 - pnbA3 /pnbR)
range(relD310 <- 1 - pnA310/pnbR)
cbind(aarg3, pnbA3, pnA310, pnbR,
      relD2 = signif(relD2, 3), relD310 = signif(relD310, 3)) # <------> Table 3
stopifnot(abs(relD2  ) < 0.09) # max is 0.06337
stopifnot(abs(relD310) < 1e-4) # max is 3.898e-5

}
\keyword{math}
\keyword{distribution}
