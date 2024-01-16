\name{pow}
\title{X to Power of Y -- R C API \code{R_pow()}}
\alias{pow}
\alias{pow_di}
\alias{.pow}
\description{
  \code{pow(x,y)} calls \R C API \file{Rmathlib}'s \code{R_pow(x,y)}
  function to compute \code{x^y} \emph{or} when \code{try.int.y} is true
  (as by default), and \code{y} is integer valued and fits into integer
  range, \code{R_pow_di(x,y)}.

  \code{pow_di(x,y)} with integer \code{y} calls \R mathlib's \code{R_pow_di(x,y)}.
}
\usage{
pow   (x, y, try.int.y = TRUE)
pow_di(x, y)
.pow  (x, y)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{y}{a numeric or in the case of \code{pow_di()} integer vector.}
  \item{try.int.y}{logical indicating if \code{pow()} should check if
    \code{y} is integer valued and fits into integer range, and in that
    case call \code{pow_di()} automatically.}
}
\details{
  In January 2024, I found (e.g., in \file{tests/pow-tst.R}) % ../tests/pow-tst.R
  that the accuracy of \code{pow_di()}, i.e., also the C function
  \code{R_pow_di()} in \R's API is of much lower precision than \R's
  \code{x^y} or (equivalently) \code{R_pow(x,y)} in \R's API, notably on
  Linux and macOS, using glib etc, sometimes as soon as \eqn{y \ge 6}{y >= 6}
  or so.

  \code{.pow(x,y)} is identical to \code{pow(x,y, try.int.y = FALSE)}
}
\value{
  a numeric vector like \code{x} or \code{y} which are recycled to common
  length, of course.
}
%% \references{

%% }
\author{Martin Maechler}
\seealso{
  Base \R's \code{\link{^}} \dQuote{operator}.
}
\examples{% -- see also ../tests/pow-tst.R
set.seed(27)
x <- rnorm(100)
y <- 0:9
stopifnot(exprs = {
    all.equal(x^y, pow(x,y))
    all.equal(x^y, pow(x,y, FALSE))
    all.equal(x^y, pow_di(x,y))
})
}
\keyword{arith}
\keyword{math}
