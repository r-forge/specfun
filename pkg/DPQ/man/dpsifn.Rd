\name{dpsifn}
\alias{dpsifn}
\title{Psi Gamma Functions Workhorse from R's API}
\description{
  Log Gamma derivatives, Psi Gamma functions.  \code{dpsifn()} is an \R
  interface to the \R API function \code{R_dpsifn()}.
}
\usage{
dpsifn(x, m, deriv1 = 0L, k2 = FALSE)
}
\arguments{
  \item{x}{numeric vector.}
  \item{m}{number of derivatives to return, an integer >= 0.}
  \item{deriv1}{\dQuote{start} derivative ....}%FIXME
  \item{k2}{a \code{\link{logical}} specifying if \code{kode = 2} should be
    applied.}
}
\details{
  \code{dpsifn()} is the underlying \dQuote{workhorse} of \R's own
  \code{\link{digamma}}, \code{\link{trigamma}} and (generalized)
  \code{\link{psigamma}} functions.

  It is useful, e.g., when several derivatives of
  \eqn{\log\Gamma=}{log Gamma=}\code{\link{lgamma}} are desired.  It
  computes and returns length-\var{m} sequence
  \eqn{(-1)^{k+1} / \Gamma(k+1) \cdot \psi^{(k)}(x)}{(-1)^(k+1) / gamma(k+1) * psi(k,x)}
  for
  \eqn{k = n, n+1,\ldots, n+m-1}{k = n, n+1, ..., n+m-1}, where
  \eqn{n=}\code{deriv1}, and \eqn{\psi^{(k)}(x)}{psi(k,x)} is the k-th
  derivative of \eqn{\psi(x)}{Psi(x)}, i.e., \code{psigamma(x,k)}.  For
  more details, see the comments in \file{src/nmath/polygamma.c}.
}
\value{
  A numeric \eqn{l_x \times m}{lx * m} \code{\link{matrix}} (where
  \eqn{l_x=}{lx=}\code{length(x)}) of scaled \eqn{\psi^{(k)}(x)}{psi(k,x)}
  values.  The matrix has \code{\link{attributes}}
  \item{underflow}{of \eqn{l_x}{lx} integer counts of the number of under- and
    over-flows, in computing the corresponding i-th matrix column for \code{x[i]}.}
  \item{ierr}{length-\eqn{l_x}{lx} integer vector of error codes, where
    \code{0} is normal/successful.}
}
\references{
  See those in \code{\link{psigamma}}
}
\author{Martin Maechler (R interface); R Core et al., see
  \code{\link{digamma}}.
}
%% \note{
%% }
\seealso{
  \code{\link{digamma}}, \code{\link{trigamma}}, \code{\link{psigamma}}.
}
\examples{
x <- seq(-3.5, 6, by=1/4)
dpx <- dpsifn(x, m = if(getRversion() >= "4.2") 7 else 5)
dpx # in R <= 4.2.1, see that sometimes the 'nz' (under-over-flow count) was uninitialized !!
j <- -1L+seq_len(nrow(dpx)); (fj <- (-1)^(j+1)*gamma(j+1))
## mdpsi <- cbind(di =   digamma(x),      -dpx[1,],
## 	       tri=  trigamma(x),       dpx[2,],
## 	       tetra=psigamma(x,2),  -2*dpx[3,],
## 	       penta=psigamma(x,3),   6*dpx[4,],
## 	       hexa =psigamma(x,4), -24*dpx[5,],
## 	       hepta=psigamma(x,5), 120*dpx[6,],
## 	       octa =psigamma(x,6),-720*dpx[7,])
## cbind(x, ie=attr(dpx,"errorCode"), round(mdpsi, 4))
str(psig <- outer(x, j, psigamma))
dpsi <- t(fj * (`attributes<-`(dpx, list(dim=dim(dpx)))))
if(getRversion() >= "4.2") {
      print( all.equal(psig, dpsi, tol=0) )# -> see 1.185e-16
  stopifnot( all.equal(psig, dpsi, tol=1e-15) )
} else { # R <= 4.1.x; dpsifn(x, ..) *not* ok for x < 0
  i <- x >= 0
      print( all.equal(psig[i,], dpsi[i,], tol=0) )# -> see 1.95e-16
  stopifnot( all.equal(psig[i,], dpsi[i,], tol=1e-15) )
}
}
\keyword{math}
\concept{psi gamma functions}
