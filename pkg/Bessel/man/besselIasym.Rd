\name{besselIasym}
\alias{besselIasym}
\alias{besselI.ftrms}
\title{Asymptotic Expansion of besselI(x,nu) For Large x}
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
besselIasym  (x, nu, k.max = 10, expon.scaled = FALSE, log = FALSE)
besselI.ftrms(x, nu, K = 20)
}
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{nu}{ ~~Describe \code{nu} here~~ }
  \item{k.max, K}{ ~~Describe \code{k.max} here~~ }
  \item{expon.scaled}{ ~~Describe \code{expon.scaled} here~~ }
  \item{log}{ ~~Describe \code{log} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{
  Abramowitz, M., and Stegun, I. A. (1955, etc).
  \emph{Handbook of mathematical functions}
  (NBS AMS series 55, U.S. Dept. of Commerce).
}
\author{Martin Maechler}
\seealso{
  This package \code{\link{BesselI}}, \pkg{base} \code{\link{besselI}}, etc
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, nu, k.max = 10, expon.scaled=FALSE, log=FALSE)
{
    ## Purpose: Asymptotic expansion of Bessel I_nu(x) function   x -> oo
    ##        by Abramowitz & Stegun (9.7.1), p.377 :
    ##
    ## I_a(z) = exp(z) / sqrt(2*pi*z) * f(z,..)  where
    ##   f(z,..) = 1 - (mu-1)/ (8*z) + (mu-1)(mu-9)/(2! (8z)^2) - ...
    ##           = 1- (mu-1)/(8z)*(1- (mu-9)/(2(8z))*(1-(mu-25)/(3(8z))*..))
    ## where  mu = 4*a^2
    ## ----------------------------------------------------------------------
    ## Arguments: x, nu, expon.scaled:  as besselI()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 11 Apr, 21 Nov 2008

    ## Note that for each x the series eventually *DIVERGES*
    ## it should be "stopped" as soon as it *converged* (:-)

    stopifnot(k.max == round(k.max))

    ## I_\nu(z) = exp(z) / sqrt(2*pi*z) * f(z, nu)

    ## First compute  f(x, nu) to order k.max
    f <- 1
    d <- 0 ## d = 1 - f  <==>  f = 1 - d
    if(k.max >= 1) {
        m. <- 4*nu^2
        x8 <- 8*x
        for(k in k.max:1) {
            ## mu <- 4*nu^2; d <- (1 - d)*(mu - (2*k-1)^2)/(k*x8)
            ## mu - (2k-1)^2 = (2nu - (2k-1)) (2nu + (2k-1)) =
            ##               = (2(nu-k)+1)(2(nu+k)-1)
            d <- (1 - d)*((2*(nu-k)+1)*(2*(nu+k)-1))/(k*x8)
        }
    }
    if(log) {
        ## f = 1 - d  ==> log(f) = log1p(-d) :
        (if(expon.scaled) log1p(-d) else x + log1p(-d)) - log(2*pi*x) / 2
    } else {
        (if(expon.scaled) (1-d) else exp(x)*(1-d)) / sqrt(2*pi*x)
    }
  }
}
\keyword{math}
