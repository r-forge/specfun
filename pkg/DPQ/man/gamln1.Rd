\name{gamln1}
\title{Compute  log( Gamma(x+1) ) Accurately in [-0.2, 1.25]}
\alias{gamln1.}
\alias{gamln1}
\description{
  Computes \eqn{\log \Gamma(a+1)} accurately notably when
  \eqn{|a| \ll 1}{|a| << 1}.
  Specifically, in the \code{userTOMS} case, it uses high (double
  precision) accuracy rational approximations for
  \eqn{-0.2 \le a \le 1.25}.
}
\usage{
gamln1 (a, useTOMS = is.numeric(a), warnIf = TRUE)
gamln1.(a, warnIf = TRUE)
}
\arguments{
  \item{a}{a numeric or numeric-alike, typically inheriting from \code{class} \code{"mpfr"}.}
  \item{useTOMS}{\code{logical} indicating if the TOMS 708 C code should be
    used; in that case \code{a} must be (C-level coercible to) numeric.}
  \item{warnIf}{logical if a \code{\link{warning}} should be signalled when
    \code{a} is not in the \dQuote{proper} range \eqn{[-0.2, 1.25]}.}
  %% \item{useDbl}{logical determining if the double precision approximation
  %%   in \eqn{[-0.2, 1.25]} from TOMS 708 should be used.}
}
\details{
  It uses \eqn{-a * p(a)/q(a)} for \eqn{a < 0.6}, where \eqn{p} and \eqn{q} are
  polynomials of degree 6 with coefficient vectors \eqn{p = [p_0 p_1 \dots p_6]}
  and \eqn{q}, \preformatted{
    p <- c( .577215664901533, .844203922187225, -.168860593646662,
	    -.780427615533591, -.402055799310489, -.0673562214325671,
	    -.00271935708322958)
    q <- c( 1, 2.88743195473681, 3.12755088914843, 1.56875193295039,
	      .361951990101499, .0325038868253937, 6.67465618796164e-4)
  }

  Similarly, for \eqn{a \ge 0.6}, \eqn{x := a - 1}, the result is
  \eqn{x * r(x)/s(x)}, with 5th degree polynomials \eqn{r()} and \eqn{s()}
  and coefficient vectors \preformatted{
    r <- c(.422784335098467, .848044614534529, .565221050691933,
           .156513060486551, .017050248402265, 4.97958207639485e-4)
    s <- c( 1 , 1.24313399877507, .548042109832463,
           .10155218743983, .00713309612391, 1.16165475989616e-4)
  }
}
\value{
  a numeric-alike vector like \code{a}.
}
\references{
  TOMS 708, see \code{\link{pbeta}}
}
\author{Martin Maechler building on C code of TOMS 708}
%% \note{

%% }
\seealso{
  \code{\link{lgamma1p}()} for different algorithms to compute \eqn{\log \Gamma(a+1)},
  notably when outside the interval \eqn{[-0.2, 1.35]}.
  Package \CRANpkg{DPQmpfr}'s \code{\link[DPQmpfr]{lgamma1pM}()} provides
  very precise such computations.
  \code{\link{lgamma}()} (and \code{gamma()} (same page)).
}
\examples{
lg1 <- function(u) lgamma(u+1) # the simple direct form
## The curve, zeros at  u=0 & u=1:
curve(lg1, -.2, 1.25, col=2, lwd=2, n=999)
title("lgamma(x + 1)"); abline(h=0, v=0:1, lty=3)

u <- (-16:100)/80 ; set.seed(1); u <- sample(u) # permuted (to check logic)
g11   <- vapply(u, gamln1, numeric(1))
gamln1. <- gamln1(u)
stopifnot( identical(g11, gamln1.) )
stopifnot( all.equal(lg1(u), g11) )

u <- (-160:1000)/800
relE <- sfsmisc::relErrV(gamln1(u), lg1(u))
plot(u, relE, type="l", main = expression("rel.diff." ~~ gamln1(u) \%~~\% lgamma(u+1)))%
plot(u, abs(relE), type="l", log="y", yaxt="n",
     main = expression("|rel.diff.|" ~~ gamln1(u) \%~~\% lgamma(u+1)))%
sfsmisc::eaxis(2)


if(requireNamespace("DPQmpfr")) withAutoprint({
  ## Comparison using Rmpfr; extending the [-.2, 1.25] interval a bit
  u <- seq(-0.225, 1.31, length.out = 2000)
  lg1pM <- DPQmpfr::lgamma1pM(Rmpfr::mpfr(u, 128))
  relE <- Rmpfr::asNumeric(sfsmisc::relErrV(lg1pM, gamln1(u, warnIf=FALSE)))

  plot(relE ~ u, type="l", ylim = c(-1,1) * 2.3e-15,
       main = expression("relative error of " ~~ gamln1(u) == log( Gamma(u+1) )))
  grid(lty = 3); abline(v = c(-.2, 1.25), col = adjustcolor(4, 1/2), lty=2, lwd=2)
  ## well... TOMS 708 gamln1() is good (if "only" 14 digits required

  ## what about the direct formula -- how bad is it really ?
  relED <- Rmpfr::asNumeric(sfsmisc::relErrV(lg1pM, lg1(u)))
  lines(relED ~ u, col = adjustcolor(2, 1/2))
  ## amazingly, the direct formula is partly (around -0.2 and +0.4) even better than gamln1() !

  plot(abs(relE) ~ u, type="l", log = "y", ylim = c(7e-17, 1e-14),
       main = expression("|relative error| of " ~~ gamln1(u) == log( Gamma(u+1) )))
  grid(lty = 3); abline(v = c(-.2, 1.25), col = adjustcolor(4, 1/2), lty=2, lwd=2)
  relED <- Rmpfr::asNumeric(sfsmisc::relErrV(lg1pM, lg1(u)))
  lines(abs(relED) ~ u, col = adjustcolor(2, 1/2))
})
}
\keyword{math}
