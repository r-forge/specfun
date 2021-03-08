\name{pnormAsymp}
\alias{pnormAsymp}
\title{Asymptotic Approxmation of (Extreme Tail) 'pnorm()'}
\description{
  Provide the first few terms of the asymptotic series approximation to
  \code{\link{pnorm}()}'s (extreme) tail, from Abramawitz and Stegun's
  26.2.13 (p.932).% \url{https://people.math.sfu.ca/~cbm/aands/page_932.htm}.
}
\usage{
pnormAsymp(x, k, lower.tail = FALSE, log.p = FALSE)
}
\arguments{
  \item{x}{positive (at least non-negative) numeric vector.}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{pnorm}()}.}
  \item{k}{integer \eqn{\ge 0}{>= 0} indicating how many terms the approximation
    should use; currently \eqn{k \le 5}{k <= 5}.}
}
\value{
  a numeric vector \dQuote{as} \code{x}; see the examples, on how to use it
  with arbitrary precise \code{\link[Rmpfr]{mpfr}}-numbers from package \CRANpkg{Rmpfr}.
}
\references{
  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain.
}
\author{Martin Maechler}
\seealso{
  \code{\link{pnormU_S53}} for (also asymptotic) upper and lower bounds.
}
\examples{
x <- c((2:10)*2, 25, (3:9)*10, (1:9)*100, (1:8)*1000, (2:4)*5000)
Px <- pnorm(x, lower.tail = FALSE, log.p=TRUE)
PxA <- sapply(setNames(0:5, paste("k =",0:5)),
              pnormAsymp, x=x, lower.tail = FALSE, log.p=TRUE)
## rel.errors :
signif(head( cbind(x, 1 - PxA/Px) , 20))

## Look more closely with high precision computations
if(requireNamespace("Rmpfr")) {
  ## ensure our function uses Rmpfr's dnorm(), etc:
  environment(pnormAsymp) <- asNamespace("Rmpfr")
  environment(pnormU_S53) <- asNamespace("Rmpfr")
  x. <- Rmpfr::mpfr(x, precBits=256)
  Px. <- Rmpfr::pnorm(x., lower.tail = FALSE, log.p=TRUE)
  ## manual, better sapplyMpfr():
  PxA. <- sapply(setNames(0:5, paste("k =",0:5)),
                 pnormAsymp, x=x., lower.tail = FALSE, log.p=TRUE)
  PxA. <- new("mpfrMatrix", unlist(PxA.), Dim=dim(PxA.), Dimnames=dimnames(PxA.))
  PxA2 <- Rmpfr::cbind(pn_dbl = Px, PxA.,
                       pnormU53 = pnormU_S53(x=x., lower.tail = FALSE, log.p=TRUE))
  ## rel.errors :
  print( Rmpfr::roundMpfr(Rmpfr::cbind(x., 1 - PxA2/Px.), precBits = 13) )
  pch <- c("R", 0:5, "U")
  matplot(x, abs(1 -PxA2/Px.), type="o", log="xy", pch=pch,
          main="pnorm(<tail>) approximations' relative errors")
  legend("bottomleft", colnames(PxA2), col=1:6, pch=pch, lty=1:5, bty="n", inset=.01)
  at1 <- axTicks(1, axp=c(par("xaxp")[1:2], 3))
  axis(1, at=at1)
  abline(h = 1:2* 2^-53, v = at1, lty=3, col=adjustcolor("gray20", 1/2))
  axis(4, las=2, at= 2^-53, label = quote(epsilon[C]), col="gray20")
}
}
\keyword{distribution}
