\name{gamln1}
\alias{gamln1}
\title{Compute  log( Gamma(x+1) ) Accurately in [-0.2, 1.25]}
\description{
  Computes \eqn{\log \Gamma(a+1)} accurately notably when
  \eqn{|a| \ll 1}{|a| << 1}.
  Specifically, it uses high (double precision) accuracy rational
  approximations for  \eqn{-0.2 \le a \le 1.25}.
}
\usage{
gamln1(a, useDbl = is.numeric(a))
}
\arguments{
  \item{a}{a numeric or numeric-alike, typically inheriting from \code{class} \code{"mpfr"}.}
  \item{useDbl}{logical determining if the double precision approximation
    in \eqn{[-0.5, 1.5]} from TOMS 708 should be used.}
  %% \item{verbose}{ ..}
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
  TODO
}
\author{Martin Maechler building on C code of TOMS 708}
%% \note{

%% }
\seealso{
  \code{\link{lgamma1p}()} for different algorithms to compute \eqn{\log \Gamma(a+1)};
  %% FIXME (in future): DPQmpfr's  lgamma1pM()
  \code{\link{lgamma}()} (and \code{gamma()} (same page)).
}
\examples{
lg1 <- function(u) lgamma(u+1) # the simple direct form
u <- seq(-.5, 1.5, by=1/16); set.seed(1); u <- sample(u) # permuted (to check logic)
%% u <- seq(-.5, 1.5, by=1/4) ## for testing
g11   <- vapply(u, gamln1, numeric(1))
gamln1. <- gamln1(u)
stopifnot( all.equal(lg1(u), g11) )
stopifnot( identical(g11, gamln1.) )

## Comparison using Rmpfr; slightly extending the [-.5, 1.5] interval:
u <- seq(-0.525, 1.525, length.out = 2001)
lg1pM <- lgamma1p(mpfr(u, 128))
relE <- asNumeric(relErrV(lg1pM, gamln1(u)))

plot(relE ~ u, type="l", ylim = c(-1,1) * 2.5e-15,
     main = expression("relative error of " ~~ gamln1(u) == log( Gamma(u+1) )))
grid(lty = 3); abline(v = c(-.5, 1.5), col = adjustcolor(4, 1/2), lty=2, lwd=2)

## what about the direct formula -- how bad is it really ?
relED <- asNumeric(relErrV(lg1pM, lg1(u)))
lines(relED ~ u, col = adjustcolor(2, 1/2), lwd = 2)
}
\keyword{math}
