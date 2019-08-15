\name{qgammaAppr}
\title{Compute (Approximate) Quantiles of the Gamma Distribution}
\alias{qgammaAppr}
\alias{qgamma.R}
\alias{qgammaApprKG}
\alias{qgammaApprSmallP}
\alias{.qgammaApprBnd}
\description{
  Compute approximations to the quantile (i.e., inverse cumulative)
  function of the Gamma distribution.
}
\usage{
qgammaAppr(p, shape, lower.tail = TRUE, log.p = FALSE,
           tol = 5e-07)
qgamma.R  (p, alpha, scale = 1, lower.tail = TRUE, log.p = FALSE,
           EPS1 = 0.01, EPS2 = 5e-07, epsN = 1e-15, maxit = 1000,
           pMin = 1e-100, pMax = (1 - 1e-14),
           verbose = getOption("verbose"))

qgammaApprKG(p, shape, lower.tail = TRUE, log.p = FALSE)
 %% simply uses 0.5* qchisqKG(p, df=2*shape, ..)

qgammaApprSmallP(p, shape, lower.tail = TRUE, log.p = FALSE)
 %% qgamma() approximation for small p -- particularly useful for small shape !
 %% based on \code{\link{lgamma1p}()}

}
\arguments{
  \item{p}{numeric vector (possibly log tranformed) probabilities.}
  \item{shape, alpha}{shape parameter, non-negative.}
  \item{scale}{scale parameter, non-negative, see \code{\link{qgamma}}.}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qgamma}()}; must
    have length 1.}

  \item{tol}{tolerance of maximal approximation error.}
  \item{EPS1}{small positive number. ...}% TODO
  \item{EPS2}{small positive number. ...}% TODO
  \item{epsN}{small positive number. ...}% TODO
  \item{maxit}{maximal number of iterations. ...}% TODO
  \item{pMin, pMax}{boundaries for \code{p}. ...}% TODO
  \item{verbose}{logical indicating if the algorithm should produce
    \dQuote{monitoring} information.}
}
\details{

  \code{qgammaApprSmallP(p, a)} should be a good approximation in the
  following situation when both \code{p} and \code{shape} \eqn{= \alpha =:
  a} are small :

 If we look at  Abramowitz&Stegun \eqn{gamma*(a,x) = x^-a * P(a,x)}
 and its series  \eqn{g*(a,x) = 1/gamma(a) * (1/a - 1/(a+1) * x + ...)},

 then the first order approximation \eqn{P(a,x) = x^a * g*(a,x) ~= x^a/gamma(a+1)}
 and hence its inverse  \eqn{x = qgamma(p, a) ~= (p * gamma(a+1)) ^ (1/a)}
 should be good as soon as  \eqn{1/a >> 1/(a+1) * x}

 <==>  x << (a+1)/a = (1 + 1/a)

 <==>  x < eps *(a+1)/a

 <==>  log(x) < log(eps) + log( (a+1)/a ) = log(eps) + log((a+1)/a)  ~  -36 - log(a)
     where log(x) ~= log(p * gamma(a+1)) / a = (log(p) + lgamma1p(a))/a

 such that the above

 <==>  (log(p) + lgamma1p(a))/a < log(eps) + log((a+1)/a)

 <==>  log(p) + lgamma1p(a) < a*(-log(a)+ log(eps) + log1p(a))

 <==>  log(p) <  a*(-log(a)+ log(eps) + log1p(a)) - lgamma1p(a) =: bnd(a)

 Note that \code{qgammaApprSmallP()} indeed also builds on \code{\link{lgamma1p}()}.

 \code{.qgammaApprBnd(a)} provides this bound \eqn{bnd(a)};
 it is simply \code{a*(logEps + log1p(a) - log(a)) - lgamma1p(a)}, where
 \code{logEps} is \eqn{\log(\epsilon)}{log(eps)} = \code{log(eps)} where \code{eps <-
 .Machine$double.eps}, i.e. typically (always?) \code{logEps}\eqn{= \log
 \epsilon = -52 * \log(2) =  -36.04365}.
}

\value{
  numeric
}
\references{
  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain.\cr
  %% Formula (???), p.???
}
\author{Martin Maechler}
\seealso{
  \code{\link{qgamma}} for \R's Gamma distribution functions.
}
\examples{
  ## TODO :  Move some of the curve()s from ../tests/qgamma-ex.R !!
}
\keyword{distribution}
