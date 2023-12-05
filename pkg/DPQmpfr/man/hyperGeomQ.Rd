\name{dhyperQ}
\title{Exact Hypergeometric Distribution Probabilites}
\alias{dhyperQ}
\alias{phyperQ}
\alias{phyperQall}
\description{% -- ../R/hyperDist.R
  Computes \bold{exact} probabilities for the hypergeometric distribution
  (see, e.g., \code{\link{dhyper}()} in \R), using package \CRANpkg{gmp}'s
  big integer and rational numbers, notably \code{\link[gmp]{chooseZ}()}.
}
\usage{
dhyperQ(x, m, n, k)
phyperQ(x, m, n, k, lower.tail=TRUE)
phyperQall(m, n, k, lower.tail=TRUE)
}
\arguments{%% >> R's ~/R/D/r-devel/R/src/library/stats/man/Hypergeometric.Rd
  \item{x}{the number of white balls drawn without replacement from an urn
    which contains both black and white balls.}
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in
    \eqn{0,1,\dots, m+n}.}
  \item{lower.tail}{logical indicating if the lower or upper tail
    probability should be computed.}
}
%% \details{
%% }
\value{
  a bigrational (class \code{"bigq"} from package \CRANpkg{gmp}) vector
  \dQuote{as} \code{x}; currently of length one (as all the function
  arguments must be \dQuote{scalar}, currently).
}
\author{Martin Maechler}
\seealso{
  \code{\link[gmp]{chooseZ}} (pkg \CRANpkg{gmp}),
  and \R's own \code{\link[stats]{Hypergeometric}}
}
\examples{
## dhyperQ() is simply
 function (x, m, n, k)
 {
    stopifnot(k - x == as.integer(k - x))
    chooseZ(m, x) * chooseZ(n, k - x) / chooseZ(m + n, k)
 }

# a case where  phyper(11, 15, 0, 12, log=TRUE) gave 'NaN'
(phyp5.0.12 <- cumsum(dhyperQ(0:12, m=15,n=0,k=12)))
stopifnot(phyp5.0.12 == c(rep(0, 12), 1))

for(x in 0:9)
  stopifnot(phyperQ(x, 10,7,8) +
            phyperQ(x, 10,7,8, lower.tail=FALSE) == 1)

(ph. <- phyperQall(m=10, n=7, k=8))
## Big Rational ('bigq') object of length 8:
## [1] 1/2431    5/374     569/4862  2039/4862 3803/4862 4685/4862 4853/4862 1
stopifnot(identical(gmp::c_bigq(list(0, ph.)),
                    1- c(phyperQall(10,7,8, lower.tail=FALSE), 0)))

(doExtras <- DPQmpfr:::doExtras())
if(doExtras) { # too slow for standard testing
 k <- 5000
 system.time(ph <-   phyper(k, 2*k, 2*k, 2*k)) #   0 (< 0.001 sec)
 system.time(phQ <- phyperQ(k, 2*k, 2*k, 2*k)) # 5.6 (was 6.3) sec
 ## Relative error of R's phyper()
 stopifnot(print(gmp::asNumeric(1 - ph/phQ)) < 1e-14) # seen 1.063e-15
}
}
\keyword{ distribution }
