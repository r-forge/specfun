\name{rexpm1}
\alias{rexpm1}
\title{TOMS 708 Approximation REXP(x) of expm1(x) = exp(x) - 1}
\description{
  Originally \code{REXP()}, now \code{rexpm1()} is a numeric (double
  precision) approximation of \eqn{exp(x) - 1},
  notably for small \eqn{|x| \ll 1}{|x| << 1} where direct evaluation
  looses accuracy through cancellation.

  Fully accurate computations of \eqn{exp(x) - 1} are now known as
  \code{\link{expm1}(x)} and have been provided by math libraries (for C,
  C++, ..) and \R, (and are typically more accurate than \code{rexp1()}).

  The \code{rexpm1()} approximation
  was developed by Didonato & Morris (1986) and uses a minimax rational
  approximation for \eqn{|x| <= 0.15}; the authors say
  \dQuote{\emph{accurate to within 2 units of the 14th significant digit}}
  (top of p.379).
}
\usage{
rexpm1(x)
}
\arguments{
  \item{x}{a numeric vector.}
}
\value{
  a numeric vector (or array) as \code{x},
}
\references{
  Didonato, A.R. and Morris, A.H. (1986)
  Computation of the Incomplete Gamma Function Ratios and their Inverse.
  \emph{ACM Trans. on Math. Softw.} \bold{12}, 377--393, \doi{10.1145/22721.23109};
  The above is the \dQuote{flesh} of \sQute{TOMS 654}:

  Didonato, A.R. and Morris, A.H. (1987)
  Algorithm 654: FORTRAN subroutines for Compute the Incomplete Gamma
  Function Ratios and their Inverse.
  \emph{ACM Transactions on Mathematical Software} \bold{13}, 318--319, \doi{10.1145/29380.214348}.
}
\author{Martin Maechler, for the C to R *vectorized* translation.
}
\seealso{
  \code{\link{pbeta}}, where the C version of \code{rexpm1()} has been used in
  several places, notably in the original TOMS 708 algorithm.
}
\examples{
x <- seq(-3/4, 3/4, by=1/1024)
plot(x,     rexpm1(x)/expm1(x) - 1, type="l", main = "Error wrt expm1()")
abline(h = (-8:8)*2^-53, lty=1:2, col=adjustcolor("gray", 1/2))
cb2 <- adjustcolor("blue", 1/2)
do.15 <- function(col = cb2) {
    abline(v = 0.15*(-1:1), lty=3, lwd=c(3,1,3), col=col)
    axis(1, at=c(-.15, .15), col=cb2, col.axis=cb2)
}
do.15()

 op <- par(mar = par("mar") + c(0,0,0,2))
plot(x, abs(rexpm1(x)/expm1(x) - 1),type="l", log = 'y',
     main = "*Relative* Error wrt expm1() [log scale]")#, yaxt="n"
abline(h = (1:9)*2^-53, lty=2, col=adjustcolor("gray", 1/2))
axis(4, at = (1:9)*2^-53, las = 1, labels =
     expression(2^-53, 2^-52, 3 \%*\% 2^-53, 2^-51, 5 \%*\% 2^-53,
                6 \%*\% 2^-53, 7 \%*\% 2^-53, 2^-50, 9 \%*\% 2^-53))
do.15()
 par(op)

## "True" Accuracy comparison of  rexpm1() with [OS mathlib based] expm1():
if(require("Rmpfr")) withAutoprint({
  xM <- mpfr(x, 128); Xexpm1 <- expm1(xM)
  REr1 <- asNumeric(rexpm1(x)/Xexpm1 - 1)
  REe1 <- asNumeric(expm1(x) /Xexpm1 - 1)
  absC <- function(E) pmax(2^-55, abs(E))

  plot(x, absC(REr1), type= "l", log="y",
       main = "|rel.Error|  of exp(x)-1 computations wrt 128-bit MPFR ")
  lines(x, absC(REe1), col = (c2 <- adjustcolor(2, 3/4)))
  abline(h = (1:9)*2^-53, lty=2, col=adjustcolor("gray60", 1/2))
  do.15()
  axis(4, mgp=c(2,1/4,0),tcl=-1/8, at=2^-(53:51), labels=expression(2^-53, 2^-52, 2^-51), las=1)
  legend("topleft", c("rexpm1(x)", " expm1(x)"), lwd=2, col=c("black", c2),
         bg = "gray90", box.lwd=.1)

})
}
\keyword{math}
\keyword{distribution}
