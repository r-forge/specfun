\name{besselI.nuAsym}
\alias{besselI.nuAsym}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
besselI.nuAsym(x, nu, k.max, expon.scaled = FALSE, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{nu}{ ~~Describe \code{nu} here~~ }
  \item{k.max}{ ~~Describe \code{k.max} here~~ }
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
function(x, nu, k.max, expon.scaled=FALSE, log=FALSE)
{
    ## Purpose: Asymptotic expansion of Bessel I_nu(x) function
    ##        when BOTH  nu and x  are large

    ## Abramowitz & Stegun , p.378, __ 9.7.7. __
    ##
    ## I_nu(nu * z) ~ 1/sqrt(2*pi*nu) * exp(nu*eta)/(1+z^2)^(1/4) *
    ##                * {1 + u_1(t)/nu + u_2(t)/nu^2 + ... }

    ## where   t = 1 / sqrt(1 + z^2),
    ##       eta = sqrt(1 + z^2) + log(z / (1 + sqrt(1+z^2)))
    ##
    ##
    ## and u_k(t)  from  p.366  __ 9.3.9 __

    ## u0(t) = 1
    ## u1(t) = (3*t - 5*t^3)/24
    ## u2(t) = (81*t^2 - 462*t^4 + 385*t^6)/1152
    ## ...

    ## with recursion  9.3.10    for  k = 0, 1, .... :
    ##
    ## u_{k+1}(t) = t^2/2 * (1 - t^2) * u'_k(t) +
    ##            1/8  \int_0^1 (1 - 5*t^2)* u_k(t) dt

    ## ----------------------------------------------------------------------
    ## Arguments: x, nu, expon.scaled:  as besselI()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 22 Nov 2008, 15:22

    stopifnot(k.max == round(k.max),
              0 <= k.max, k.max <= 3)

    z <- x/nu # -> computing I_nu(z * nu)
    sz <- sqrt(1 + z^2)
    t <- 1/sz
    eta <- (if(expon.scaled) ## <==> * exp(-x)  ==  exp(- z * nu)  <==> "eta - z"
            ## sz - z = sqrt(1 + z^2) - z  =  1/(sqrt(..) + z)
            1/(z + sz) else sc) +
                + log(z / (1 + sz))

    ## I_nu(nu * z) ~ 1/sqrt(2*pi*nu) * exp(nu*eta)/(1+z^2)^(1/4) *
    ##                * {1 + u_1(t)/nu + u_2(t)/nu^2 + ... }
    if(k.max == 0)
        d <- 0
    else { ## k.max >= 1
        t2 <- t^2
        u1.t <- (t*(3 - 5*t2))/24
        d <-
            if(k.max == 1) {
                u1.t/nu
            } else { ## k.max >= 2
                u2.t <- ##(81*t^2 - 462*t^4 + 385*t^6)/1152
                    t2*(81 + t2*(-462 + t2*385)) / 1152
                if(k.max == 2)
                    (u1.t + u2.t/nu)/nu
                else { ## k.max == 3
                    u3.t <- t*t2*(30375 +
                                  t2*(-369603 +
                                      t2*(765765 - t2*425425)))/ 414720
                    ## if(k.max == 3)
                    (u1.t + (u2.t + u3.t/nu)/nu)/nu
                }
            }
    }

    if(log) {
        log1p(d) + nu*eta - (log(sz) + log(2*pi*nu))/2
    } else {
        (1+d) * exp(nu*eta)/sqrt(2*pi*nu*sz)
    }
  }
}
\keyword{math}

