\name{pl2curves}
\alias{pl2curves}
\title{Plot 2 Noncentral Distribution Curves for Visual Comparison}
\description{
  Plot two noncentral (chi-squared or \eqn{t} or ..) distribution curves
  for visual comparison.
}
\usage{% >>>>> ../R/dnchisq-fn.R
pl2curves(fun1, fun2, df, ncp, log = FALSE,
          from = 0, to = 2 * ncp, p.log = "", n = 2001,
          leg = TRUE, col2 = 2, lwd2 = 2, lty2 = 3, \dots)
}
\arguments{
  \item{fun1, fun2}{\code{\link{function}()}s, both to be used via
    \code{\link{curve}()}, and called with the same 4 arguments,
    \code{(., df, ncp, log)} (the name of the first argument is not specified).}
  \item{df, ncp, log}{parameters to be passed and used in both functions,
    which hence typically are non-central chi-squared or t density,
    probability or quantile functions.}
  %% -- parameters for curve() :
  \item{from, to}{numbers determining the x-range, passed to
    \code{\link{curve}()}.}
  \item{p.log}{string, passed as \code{\link{curve}(...., log = log.p)}.}
  \item{n}{the number of evaluation points, passed to \code{\link{curve}()}.}
  \item{leg}{logical specifying if a \code{\link{legend}()} should be drawn.}
  \item{col2, lwd2, lty2}{color, line width and line type for the second
    curve. (The first curve uses defaults for these graphical properties.)}
  \item{\dots}{further arguments passed to \emph{first} \code{\link{curve}(..)} call.}
}
%% \details{
%% }
\value{
  TODO: inivisible return both curve() results, i.e., (x,y1, y2), possibly
  as data frame
}
\author{Martin Maechler}
\seealso{
  \code{\link{curve}}, ..
}
\examples{

p.dnchiBessel <- function(df, ncp, log=FALSE, from=0, to = 2*ncp, p.log="", ...)
{
    pl2curves(dnchisqBessel, dchisq, df=df, ncp=ncp, log=log,
              from=from, to=to, p.log=p.log, ...)
}

  ## TODO the p.dnchiB()  examples  >>>>>> ../tests/chisq-nonc-ex.R <<<
}
\keyword{hplot}
