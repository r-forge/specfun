\name{logcf}
\alias{logcf}
\title{Continued Fraction Approximation of Log-Related Series}
\description{
  Compute a continued fraction approximation to the series (infinite sum)
  \deqn{\sum_{k=0}^\infty \frac{x^k}{i +k\cdot d} = \frac{1}{i} + \frac{x}{i+d} +
    \frac{x^2}{i+2*d} + \frac{x^3}{i+3*d} + \ldots}{%
    sum(k=0,...,Inf; x^k/(i+k*d)) = 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ...}

   Needed as auxiliary function in \code{\link{log1pmx}()} and \code{\link{lgamma1p}()}.
}
\usage{
logcf(x, i, d, eps, maxit = 10000)
}
\arguments{
  \item{x}{numeric vector}
  \item{i}{positive numeric}
  \item{d}{non-negative numeric}
  \item{eps}{positive number, the convergence tolerance.}
  \item{maxit}{a positive integer, the maximal number of iterations or
    terms in the truncated series used.}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  a numeric vector with the same attributes as \code{x}.
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{Martin Maechler, based on .... }
\note{
 Rescaling is done by (namespace hidden) \dQuote{global}
 \code{scalefactor} ....% FIXME
}
\seealso{
  \code{\link{lgamma1p}}, \code{\link{log1pmx}}, and
  \code{\link{pbeta}}, whose prinicipal algorithm has evolved from TOMS 708.
}
\examples{
l32 <- curve(logcf(x, 3,2, eps=1e-7), -3, 1)
abline(h=0,v=1, lty=3, col="gray50")
plot(y~x, l32, log="y", type = "o", main = "logcf(*, 3,2)  in log-scale")
}
\keyword{math}% one of  RShowDoc("KEYWORDS")
