\name{pow1p}
\alias{pow1p}
\title{Accurate \eqn{(1+x)^y}, notably for small \eqn{|x|}}
\description{
  Compute \eqn{(1+x)^y} accurately, notably also for small \eqn{|x|}, where
  the naive formula suffers from cancellation, returning \code{1}, often.
}
\usage{
pow1p(x, y,
      pow = ((x + 1) - 1) == x || abs(x) > 0.5 || is.na(x))
}
\arguments{
  \item{x, y}{numeric or number-like; in the latter case, arithmetic incl. \code{^},
    comparison, \code{\link{exp}}, \code{\link{log1p}}, \code{\link{abs}},
    and \code{\link{is.na}} methods must work.}
  \item{pow}{\code{\link{logical}} indicating if the \dQuote{naive} /
    direct computation \code{(1 + x)^y} should be used (unless \code{y} is
    in 0:4, where the binomial is used, see \sQuote{Details}).
    The current default is the one used in R's C-level function (but beware
    of compiler optimization there!).}
}
\details{
  A pure \R-implementation of R 4.4.0's new C-level
  \code{pow1p()} function which was introduced for more accurate
  \code{\link{dbinom_raw}()} computations.

  Currently, we use the \dQuote{exact} (nested) polynomial formula for
  \eqn{y \in \{0,1,2,3,4\}}{y = 0, 1, 2, 3, 4}.

  MM is conjecturing that the default \code{pow=FALSE} for (most)
  \eqn{x \le \frac 1 2}{x <= 1/2} is sub-optimal.
}
\value{
  numeric or number-like, as \code{x + y}.
}
%% \references{

%% }
\author{
  Originally proposed by Morten Welinder, see \PR{18642};
  tweaked, notably for small integer \code{y}, by Martin Maechler.
}
\seealso{
  \code{\link{^}}, \code{\link{log1p}},
  \code{\link{dbinom_raw}}.
}
\examples{
x <- 2^-(1:50)
y <- 99
f1 <- (1+x)^99
f2 <- exp(y * log1p(x))
fp <- pow1p(x, 99)
matplot(x, cbind(f1, f2, fp), type = "l", col = 2:4)
legend("top", legend = expression((1+x)^99, exp(99 * log1p(x)), pow1p(x, 99)),
       bty="n", col=2:4, lwd=2)
cbind(x, f1, f2, sfsmisc::relErrV(f2, f1))
}
\keyword{arith}
\keyword{math}
