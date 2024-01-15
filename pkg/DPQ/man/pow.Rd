\name{pow}
\alias{pow}
\alias{pow_di}
\title{X to Power of Y -- R C API \code{R_pow()}}
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
}
\arguments{
  \item{x}{a numeric vector.}
  \item{y}{a numeric or in the case of \code{pow_di()} integer vector.}
  \item{try.int.y}{logical indicating if \code{pow()} should check if
    \code{y} is integer valued and fits into integer range, and in that
    case call \code{pow_di()} automatically.}
}
%% \details{
%% }
\value{
  a numeric vector like \code{x} or \code{y}.
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
