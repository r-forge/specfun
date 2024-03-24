\name{dgamma.R}
\alias{dgamma.R}
\title{Gamma Density Function Alternatives}
\description{
  \code{dgamma.R()} is aimed to be an R level \dQuote{clone} of \R's C
  level implementation \code{\link{dgamma}} (from package \pkg{stats}).
}
\usage{
dgamma.R(x, shape, scale = 1, log,
         dpois_r_args = list())
}
\arguments{
  \item{x}{non-negative numeric vector.}
  \item{shape}{non-negative shape parameter of the Gamma distribution.}
  \item{scale}{positive scale parameter; note we do not see the need to have
    a \code{rate} parameter as the standard \R function.}
  \item{log}{logical indicating if the result is desired on the log scale.}
  \item{dpois_r_args}{a \code{\link{list}} of optional arguments for
    \code{\link{dpois_raw}()}; not much checked, must be specified correctly.}
}
%% \details{
%% }
\value{
  numeric vector of the same length as \code{x} (which may have to be
  thought of recycled along \code{shape} and/or \code{scale}.
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{Martin Maechler}
\seealso{
  (As \R's C code) this depends crucially on the \dQuote{workhorse}
  function \code{\link{dpois_raw}()}.
}
\examples{

xy  <- curve(dgamma  (x, 12), 0,30) # R's dgamma()
xyR <- curve(dgamma.R(x, 12, dpois_r_args = list(verbose=TRUE)), add=TRUE,
             col = adjustcolor(2, 1/3), lwd=3)
stopifnot(all.equal(xy, xyR, tolerance = 4e-15)) # seen 7.12e-16
## TODO: check *vectorization* in x --> add tests/*.R				___ TODO ___


## From R's  <R>/tests/d-p-q-r-tst-2.R -- replacing dgamma() w/ dgamma.R()
## PR#17577 - dgamma(x, shape)  for shape < 1 (=> +Inf at x=0) and very small x
stopifnot(exprs = {
    all.equal(dgamma.R(2^-1027, shape = .99 , log=TRUE), 7.1127667376, tol=1e-10)
    all.equal(dgamma.R(2^-1031, shape = 1e-2, log=TRUE), 702.8889158,  tol=1e-10)
    all.equal(dgamma.R(2^-1048, shape = 1e-7, log=TRUE), 710.30007699, tol=1e-10)
    all.equal(dgamma.R(2^-1048, shape = 1e-7, scale = 1e-315, log=TRUE),
              709.96858768, tol=1e-10)
})
## R's dgamma() gave all Inf in R <= 3.6.1 [and still there in 32-bit Windows !]
}
\keyword{math}
\keyword{distribution}
