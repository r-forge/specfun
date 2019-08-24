\name{newton}
\alias{newton}
\title{Simple R level Newton Algorithm, Mostly for Didactical Reasons}
\description{
  Given the function \code{G()} and its derivative \code{g()},
  \code{newton()} uses the Newton method, starting at \code{x0},
  to find a point xp at which G is zero.  \code{G()} and \code{g()}
  may each depend on the same parameter (vector) \code{z}.

  Convergence typically happens when the stepsize becomes smaller than
  \code{eps}.

  \code{keepAll = TRUE}  to also get  the vectors of consecutive values of
  x and G(x, z);
}
\usage{%% >> ../R/qnchisq.R <<<
newton(x0, G, g, z,
       xMin = -Inf, xMax = Inf, warnRng = TRUE,
       dxMax = 1000, eps = 0.0001, maxiter = 1000L,
       warnIter = missing(maxiter) || maxiter >= 10L,
       keepAll = NA)
}
\arguments{
  \item{x0}{numeric start value.}
  \item{G, g}{must be \code{\link{function}}s, mathematically of their first
    argument, but they can accept parameters; \code{g()} must be the
    derivative of \code{G}.}
  \item{z}{parameter vector for \eqn{G()} and \eqn{g()}, to be kept fixed.}
  \item{xMin, xMax}{numbers defining the allowed range for x during the
    iterations; e.g., useful to set to \code{0} and \code{1} during quantile
    search.}
  \item{warnRng}{\code{\link{logical}} specifying if a \code{\link{warning}} should be
    signalled when start value \code{x0} is outside \code{[xMin, xMax]} and
    hence will be changed to one of the boundary values.}
  \item{dxMax}{maximal step size in \eqn{x}-space.  (The default \code{1000}
    is quite arbitrary, do set a good maximal step size yourself!)}
  \item{eps}{positive number, the \emph{absolute} convergence tolerance.}
  \item{maxiter}{positive integer, specifying the maximal number of Newton
    iterations.}
  \item{warnIter}{logical specifying if a \code{\link{warning}} should be
    signalled when the algorithm has not converged in \code{maxiter} iterations.}
  \item{keepAll}{logical specifying if the full sequence of x- and G(x,*)
    values should be kept and returned: \describe{
      \item{\code{NA},}{the default: \code{newton} returns a
	small list of final \dQuote{data}, with 4 components
	\code{x}\eqn{ = x*}, \code{G}\eqn{= G(x*, z)},
	\code{it}, and \code{converged}.}
      \item{\code{TRUE}:}{returns an extended \code{\link{list}}, in
	addition containing the vectors \code{x.vec} and \code{G.vec}.}
      \item{\code{FALSE}:}{returns only the \eqn{x*} value.}
    }
  }
}
\details{
  Because of the quadrativc convergence at the end of the Newton algorithm,
  often \eqn{x^*}{x*} satisfies approximately \eqn{ | G(x*, z)| < eps^2 }.

  \code{newton()} can be used to compute the quantile function of a
  distribution, if you have a good starting value, and provide
  the cumulative probability and density functions as \R functions \code{G}
  and \code{g} respectively.
}
\value{
  The result always contains the final x-value \eqn{x*}, and typically some
  information about convergence, depending on the value of \code{keepAll},
  see above:
  \item{x}{the optimal \eqn{x^*}{x*} value (a number).}
  \item{G}{the function value \eqn{G(x*, z)}, typically very close to zero.}
  \item{it}{the integer number of iterations used.}
  \item{convergence}{logical indicating if the Newton algorithm converged
    within \code{maxiter} iterations.}

  \item{x.vec}{the full vector of x values, \eqn{\{x0,\ldots,x^*\}}{{x0,
	.., x*}}.}
  \item{G.vec}{the vector of function values (typically tending to zero),
    i.e., \code{G(x.vec, .)} (even when \code{G(x, .)} would not vectorize).}
}
\references{
  Newton's Method on Wikipedia,
  \url{https://en.wikipedia.org/wiki/Newton\%27s_method}.
}
\author{Martin Maechler, ca. 2004}
\seealso{
  \code{\link{uniroot}()} is much more sophisticated, works without
  derivatives and is generally faster than \code{newton()}.

  \code{newton(.)} is currently crucially used (only) in our function
  \code{\link{qchisqN}()}.
}
\examples{
## The most simple non-trivial case :  Computing SQRT(a)
  G <- function(x, a) x^2 - a
  g <- function(x, a) 2*x

  newton(1, G, g, z = 4  ) # z = a -- converges immediately
  newton(1, G, g, z = 400) # bad start, needs longer to converge

## More interesting, and related to non-central (chisq, e.t.) computations:
## When is  x * log(x) < B,  i.e., the inverse function of G = x*log(x) :
xlx <- function(x, B) x*log(x) - B
dxlx <- function(x, B) log(x) + 1

Nxlx <- function(B) newton(B, G=xlx, g=dxlx, z=B, maxiter=Inf)$x
N1   <- function(B) newton(B, G=xlx, g=dxlx, z=B, maxiter = 1)$x
N2   <- function(B) newton(B, G=xlx, g=dxlx, z=B, maxiter = 2)$x

Bs <- c(outer(c(1,2,5), 10^(0:4)))
plot (Bs, vapply(Bs, Nxlx, pi), type = "l", log ="xy")
lines(Bs, vapply(Bs, N1  , pi), col = 2, lwd = 2, lty = 2)
lines(Bs, vapply(Bs, N2  , pi), col = 3, lwd = 3, lty = 3)

BL <- c(outer(c(1,2,5), 10^(0:6)))
plot (BL, vapply(BL, Nxlx, pi), type = "l", log ="xy")
lines(BL, BL, col="green2", lty=3)
lines(BL, vapply(BL, N1  , pi), col = 2, lwd = 2, lty = 2)
lines(BL, vapply(BL, N2  , pi), col = 3, lwd = 3, lty = 3)
## Better starting value from an approximate 1 step Newton:
iL1 <- function(B) 2*B / (log(B) + 1)
lines(BL, iL1(BL), lty=4, col="gray20") ## really better ==> use it as start

Nxlx <- function(B) newton(iL1(B), G=xlx, g=dxlx, z=B, maxiter=Inf)$x
N1   <- function(B) newton(iL1(B), G=xlx, g=dxlx, z=B, maxiter = 1)$x
N2   <- function(B) newton(iL1(B), G=xlx, g=dxlx, z=B, maxiter = 2)$x

plot (BL, vapply(BL, Nxlx, pi), type = "o", log ="xy")
lines(BL, iL1(BL),  lty=4, col="gray20")
lines(BL, vapply(BL, N1  , pi), type = "o", col = 2, lwd = 2, lty = 2)
lines(BL, vapply(BL, N2  , pi), type = "o", col = 3, lwd = 2, lty = 3)
## Manual 2-step Newton
iL2 <- function(B) { lB <- log(B) ; B*(lB+1) / (lB * (lB - log(lB) + 1)) }
lines(BL, iL2(BL), col = adjustcolor("sky blue", 0.6), lwd=6)
##==>  iL2() is very close to true curve
## relative error:
iLtrue <- vapply(BL, Nxlx, pi)
cbind(BL, iLtrue, iL2=iL2(BL), relErL2 = 1-iL2(BL)/iLtrue)
## absolute error (in log-log scale; always positive!):
plot(BL, iL2(BL) - iLtrue, type = "o", log="xy", axes=FALSE)
if(requireNamespace("sfsmisc")) {
  sfsmisc::eaxis(1)
  sfsmisc::eaxis(2, sub10=2)
} else {
  cat("no 'sfsmisc' package; maybe  install.packages(\"sfsmisc\")  ?\n")
  axis(1); axis(2)
}
## 1 step from iL2()  seems quite good:
B. <- BL[-1] # starts at 2
NL2 <- lapply(B., function(B) newton(iL2(B), G=xlx, g=dxlx, z=B, maxiter=1))
str(NL2)
iL3 <- sapply(NL2, `[[`, "x")
cbind(B., iLtrue[-1], iL2=iL2(B.), iL3, relE.3 = 1- iL3/iLtrue[-1])
x. <- iL2(B.)
all.equal(iL3, x. - xlx(x., B.) / dxlx(x.)) ## 7.471802e-8
## Algebraic simplification of one newton step :
all.equal((x.+B.)/(log(x.)+1), x. - xlx(x., B.) / dxlx(x.), tol = 4e-16)
iN1 <- function(x, B) (x+B) / (log(x) + 1)
B <- 12345
iN1(iN1(iN1(B, B),B),B)
Nxlx(B)
}
\keyword{math}
