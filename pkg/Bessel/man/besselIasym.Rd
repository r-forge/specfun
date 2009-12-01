\name{besselIasym}
\alias{besselIasym}
%\alias{besselKasym}
\alias{besselI.ftrms}
\title{Asymptotic Expansion of besselI(x,nu) For Large x}
\description{
  Compute Bessel function \eqn{I_{\nu}(x)}{I[nu](x)}
  and \eqn{K_{\nu}(x)}{K[nu](x)} for large \eqn{x} and small or moderate
  \eqn{\nu}{nu}, using the asymptotic expansion (9.7.1), p.377 of
  Abramowitz & Stegun, for \eqn{x \to\infty}{x -> Inf}, even valid for
  \emph{\code{\link{complex}}} \eqn{x},
  %% CARE: With the new parser you should only carefully wrap these lines:
  \deqn{I_a(x) = exp(x) / \sqrt{2\pi x} \cdot f(x, a),
  }{I_a(x) = exp(x) / sqrt(2*pi*x)* f(x, a),}
  where
  \deqn{f(x,a) = 1 - \frac{\mu-1}{8x} + \frac{(\mu-1)(\mu-9)}{2! (8x)^2}
    - \ldots,}{f(x,a) = 1 - (mu-1) / (8x) + (mu-1)(mu-9) / (2! (8x)^2) - ...,}
  %            = 1- (mu-1)/(8z)*(1- (mu-9)/(2(8z))*(1-(mu-25)/(3(8z))*..))
  and  \eqn{\mu = 4 a^2}{mu = 4*a^2} and \eqn{|arg(x)| < \pi/2}.

  Whereas \code{besselIasym(x,a)} computes \eqn{I_a(x)},
  \code{besselI.ftrms} returns the corresponding \emph{terms} in the
  series expansion of \eqn{f(x,a)} above.
}
\usage{
besselIasym  (x, nu, k.max = 10, expon.scaled = FALSE, log = FALSE)
besselI.ftrms(x, nu, K = 20)
}
\arguments{
  \item{x}{numeric, \eqn{\ge 0}{>= 0}.}
  \item{nu}{numeric; The \emph{order} (maybe fractional!) of the
    corresponding Bessel function.}
  \item{k.max, K}{integer number of terms in the expansion.}% << TODO more ?
  \item{expon.scaled}{logical; if \code{TRUE}, the results are
    exponentially scaled in order to avoid overflow.}
  % (\eqn{I_{\nu}}{I(nu)}) or underflow (\eqn{K_{\nu}}{K(nu)}),
  % respectively.}
  \item{log}{logical; if TRUE, \eqn{\log(f(.))} is returned instead of \eqn{f}.}
}
\details{
 ......... FIXME ...
}
\value{
 a numeric vector of the same length as \code{x}.%% FIXME, x = matrix ?
}
\references{
  Abramowitz, M., and Stegun, I. A. (1955, etc).
  \emph{Handbook of mathematical functions}
  (NBS AMS series 55, U.S. Dept. of Commerce).
}
\author{Martin Maechler}
\seealso{
  From this package \pkg{Bessel()} \code{\link{BesselI}()}; further,
  \code{\link{besselI.nuAsym}()} which is useful when \eqn{\nu} is large
  (as well); further \pkg{base} \code{\link{besselI}}, etc
}
\examples{
x <- c(1:10, 20, 50, 100^(2:10))
nu <- c(1, 10, 20, 50, 100)
%% FIXME: do something more sensible -- compare with  besselI() and *K()
r <- lapply(c(0:4,10,20), function(k.)
            sapply(nu, function(n.)
                  besselIasym(x, nu=n., k.max = k., log = TRUE)))
warnings() %% NaN's  ... seem a bit scattered.  FIXME ??
}
\keyword{math}
