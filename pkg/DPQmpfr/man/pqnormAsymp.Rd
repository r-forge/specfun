\name{pqnormAsymp}
\title{Asymptotic Approximations of Extreme Tail 'pnorm()' and 'qnorm()'}
\alias{pnormAsymp}
\alias{qnormAsymp}
\description{
  These functions provide the first terms of asymptotic series approximations to
  \code{\link{pnorm}()}'s (extreme) tail, from Abramawitz and Stegun's
  26.2.13 (p.932),% \url{https://personal.math.ubc.ca/~cbm/aands/page_932.htm}
  or \code{\link{qnorm}()} where the approximations have been derived via
  iterative plugin using Abramowitz and Stegun's formula.
}
\usage{% ../R/dpq-mpfrized.R  ../../DPQ/man/pnormAsymp.Rd  and  ../../DPQ/man/qnormAsymp.Rd
pnormAsymp(x, k, lower.tail = FALSE, log.p = FALSE)
qnormAsymp(p, lp = .DT_Clog(p, lower.tail = lower.tail, log.p = log.p),
           order, M_2PI =,
           lower.tail = TRUE, log.p = missing(p))
}
\arguments{
  % ../../DPQ/man/pnormAsymp.Rd :
  \item{x}{positive (at least non-negative) numeric vector.}
  \item{k}{integer \eqn{\ge 0}{>= 0} indicating how many terms the approximation
    should use; currently \eqn{k \le 5}{k <= 5}.}
  % ../../DPQ/man/qnormAsymp.Rd :
  \item{p}{numeric vector of probabilities, possibly transformed, depending
    on \code{log.p}.  Does not need to be specified, if \code{lp} is instead.}
  \item{lp}{numeric (vector) of \code{log(1-p)} values; if not specified,
    computed from \code{p}, depending on \code{lower.tail} and \code{log.p}.}
  \item{order}{an integer in \eqn{\{0,1,\dots,5\}}{{0,1,..,5}}, specifying the
    approximation order.}
  \item{M_2PI}{the number \eqn{2\pi} in the same precision as \code{p} or
    \code{lp}, i.e., \code{\link{numeric}} or of class
    \code{"\link[Rmpfr]{mpfr}"}.}
  \item{lower.tail}{logical; if true, probabilities are \eqn{P[X \le x]},
    otherwise upper tail probabilities, \eqn{P[X > x]}.}
  \item{log.p}{logical; if \code{TRUE} (default for \code{qnormAsymp} !!), probabilities
    \eqn{p} are given as \eqn{\log(p)} in argument \code{p}
    or \eqn{\log{(1-p)}} in \code{lp}.}
}
\details{
  see \emph{both} help pages \code{\link[DPQ]{pnormAsymp}} and
  \code{\link[DPQ]{qnormAsymp}} from our package \CRANpkg{DPQ}.
}
\value{
  vector/array/mpfr like first argument \code{x} or \code{p} or \code{lp}, respectively.
}
\author{Martin Maechler}
\seealso{
  \code{\link{pnorm}}.  The same functions  \dQuote{numeric-only} are in my
 \CRANpkg{DPQ} package with more extensive documentation.
}
\examples{%% see ../../DPQ/man/pnormAsymp.Rd  and  ../../DPQ/man/pnormAsymp.Rd
require("Rmpfr") # (in strong dependencies of this pkg {DPQmpfr})
x <- seq(1/64, 10, by=1/64)
xm  <- mpfr( x, 96)
"TODO"

## More extreme tails: ----------------------------------------------
##
## 1. pnormAsymp() ---------------------
lx <- c((2:10)*2, 25, (3:9)*10, (1:9)*100, (1:8)*1000, (2:7)*5000)
lxm <- mpfr(lx, 256)
Px <- pnorm(lxm, lower.tail = FALSE, log.p=TRUE)
PxA <- sapplyMpfr(setNames(0:5, paste("k =",0:5)),
                  pnormAsymp, x=lxm, lower.tail = FALSE, log.p=TRUE)
if(interactive())
  roundMpfr(PxA, 40)
# rel.errors :
relE <- asNumeric(1 - PxA/Px)
options(width = 99) -> oop # (nicely printing the matrices)
cbind(lx, relE)
matplot(lx, abs(relE), type="b", cex = 1/2, log="xy", pch=as.character(0:5),
        axes=FALSE,
        main = "|relE( <pnormAsymp(lx, k=*, lower.tail=FALSE, log.p=TRUE) )|")
sfsmisc::eaxis(1, sub10=2); sfsmisc::eaxis(2)
legend("bottom", paste("k =", 0:5), col=1:6, lty=1:5,
       pch = as.character(0:5), pt.cex=1/2, bty="n")
## NB: rel.Errors go down to  7e-59 ==> need precision of  -log2(7e-59) ~ 193.2 bits

## 2. qnormAsymp() ---------------------
QPx <- sapplyMpfr(setNames(0:5, paste("k =",0:5)),
                  function(k) qnormAsymp(Px, order=k, lower.tail = FALSE, log.p=TRUE))
(relE.q <- asNumeric(QPx/lx - 1))
         # note how consistent the signs are (!) <==> have upper/lower bounds

matplot(-asNumeric(Px), abs(relE.q), type="b", cex = 1/2, log="xy", pch=as.character(0:5),
        xlab = quote(-Px), axes=FALSE,
        main = "|relE( <qnormAsymp(Px, k=*, lower.tail=FALSE, log.p=TRUE) )|")
sfsmisc::eaxis(1, sub10=2); sfsmisc::eaxis(2)
legend("bottom", paste("k =", 0:5), col=1:6, lty=1:5,
       pch = as.character(0:5), pt.cex=1/2, bty="n")

options(oop) # {revert to previous state}
}
\keyword{distribution}
