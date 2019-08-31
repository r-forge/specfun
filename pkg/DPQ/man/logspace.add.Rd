\name{logspace.add}
\alias{logspace.add}
\alias{logspace.sub}
\title{Logspace Arithmetix -- Addition and Subtraction}
\description{
  Compute the log(arithm) of a sum (or difference) from the log of terms
  without causing overflows and without throwing away large handfuls of accuracy.

  \describe{
    \item{\code{logspace.add(lx, ly)}:=}{
        \deqn{\log (\exp (lx) + \exp (ly))}}

    \item{\code{logspace.sub(lx, ly)}:=}{
        \deqn{\log (\exp (lx) - \exp (ly))}}
    }
}
\usage{
logspace.add(lx, ly)
logspace.sub(lx, ly)
}
\arguments{
  \item{lx, ly}{numeric vectors, typically of the same
        \code{\link{length}}, but will be recycled to common length as with
        other \R arithmetic.}
}
%% \details{

%% }
\value{
  a \code{\link{numeric}} vector of the same length as \code{x+y}.
}
%% \references{
%% }
\author{Morten Welinder (for \R's \code{pgamma()}); Martin Maechler}
\note{
  This is really from \R's C source code for \code{\link{pgamma}()}, i.e.,
  \file{<R>/src/nmath/pgamma.c}

  The function definitions are very simple, \code{logspace.sub()} using \code{\link{log1mexp}()}.
}
\seealso{
  \code{\link{lsum}}, \code{\link{lssum}}; then \code{\link{pgamma}()}
}
\examples{
set.seed(12)
ly <- rnorm(100, sd= 50)
lx <- ly + abs(rnorm(100, sd=100))  # lx - ly must be positive for *.sub()
stopifnot(exprs = {
   all.equal(logspace.add(lx,ly),
             log(exp(lx) + exp(ly)), tol=1e-14)
   all.equal(logspace.sub(lx,ly),
             log(exp(lx) - exp(ly)), tol=1e-14)
})
}
\keyword{arith}
