\name{BesselH}
\Rdversion{1.1}
\alias{Hankel}
\alias{BesselH}
\title{Hankel (H-Bessel) Function (of Complex Argument)}
\description{
  Compute the Hankel functions \eqn{H(1,*)} and \eqn{H(2,*)},
  also called \sQuote{H-Bessel} function (of the third kind),
  of complex arguments.  They are defined as
  \deqn{
    H(1,\nu, z) := H_{\nu}^{(1)}(z) =  J_{\nu}(z) + i Y_{\nu}(z),}{
    H(1, ν, z)  := H¹_ν(z) = J_ν(z) + i Y_ν(z),
  }
  \deqn{
    H(2,\nu, z) := H_{\nu}^{(2)}(z) =  J_{\nu}(z) - i Y_{\nu}(z),}{
    H(2, ν, z)  := H²_ν(z) = J_ν(z) - i Y_ν(z),
  }
  where \eqn{J_{\nu}(z)}{J__ν(z)} and \eqn{Y_{\nu}(z)}{Y_ν(z)} are the
  Bessel functions of the first and second kind, see
  \code{\link{BesselJ}}, etc.
}
\usage{
BesselH(m, z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0)
}
\arguments{
  \item{m}{integer, either 1 or 2, indicating the kind of Hankel function.}
  \item{z}{complex or numeric vector of values \bold{different from 0}.}
  \item{nu}{numeric, must currently be non-negative.}
  \item{expon.scaled}{logical indicating if the result should be scaled
    by an exponential factor (typically to avoid under- or over-flow).}
  \item{nSeq}{positive integer; if \eqn{> 1}, computes the result for
    a whole \emph{sequence} of \code{nu} values of length \code{nSeq},
    see \sQuote{Details} below.}
  \item{verbose}{integer defaulting to 0, indicating the level of verbosity
    notably from C code.}
}
\details{% ../src/zbsubs.f --> zbesh()
  By default (when \code{expon.scaled} is false), the resulting sequence
  (of length \code{nSeq}) is for \eqn{m = 1,2},
  \deqn{y_j = H(m, \nu+j-1, z),}{y[j]= H(m, nu+j-1, z),}
  computed for \eqn{j=1,\dots,nSeq}{j=1,...,nSeq}.

  If \code{expon.scaled} is true, the sequence is for \eqn{m = 1,2}
  \deqn{y_j = \exp(-\tilde{m} z i) \cdot H(m, \nu+j-1, z),}{%
        y[j]= exp(-mm*z* i)* H(m, nu+j-1, z),}
  where \eqn{\tilde{m} = 3-2m}{mm = 3-2*m}
  (and \eqn{i^2 = -1}), for \eqn{j=1,\dots,nSeq}{j=1,...,nSeq}.
}
\value{
  a complex or numeric vector (or \code{\link{matrix}} if \code{nSeq > 1})
  of the same length and \code{\link{mode}} as \code{z}.
}
\references{see \code{\link{BesselI}}.
}
\author{
  Donald E. Amos, Sandia National Laboratories, wrote the original
  fortran code.
  Martin Maechler did the \R interface.
}
\seealso{
  \code{\link{BesselI}} etc; the Airy function \code{\link{Airy}}.
}
\examples{
##------------------ H(1, *) ----------------
nus <- c(1,2,5,10)
for(i in seq_along(nus))
   curve(BesselH(1, x, nu=nus[i]), -10, 10, add= i > 1, col=i, n=1000)
legend("topleft", paste("nu = ", format(nus)), col = seq_along(nus), lty=1)

## nu = 10 looks a bit  "special" ...   hmm...
curve(BesselH(1, x, nu=10), -.3, .3, col=4,
      ylim = c(-10,10), n=1000)

##------------------ H(2, *) ----------------
for(i in seq_along(nus))
   curve(BesselH(2, x, nu=nus[i]), -10, 10, add= i > 1, col=i, n=1000)
legend("bottomright", paste("nu = ", format(nus)), col = seq_along(nus), lty=1)
## the same nu = 10 behavior ..
}
\keyword{math}
