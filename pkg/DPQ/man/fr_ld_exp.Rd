\name{fr_ld_exp}
\title{Base-2 Representation and Multiplication of Numbers}
\alias{frexp}
\alias{ldexp}
\description{
  Both are \R versions of C99 (and POSIX) standard C (and C++) mathlib
  functions of the same name.

  \code{frexp(x)} computes base-2 exponent \code{e} and \dQuote{mantissa},
  or \emph{fraction} \code{r}, such that \eqn{x = r * 2^e}, where \eqn{r \in
    [0.5, 1)} (unless when \code{x} is in \code{c(0, -Inf, Inf, NaN)}
  where \code{r == x} and \code{e} is 0),
  and \eqn{e} is integer valued.

  \code{ldexp(f, E)} is the \emph{inverse} of \code{frexp()}: Given
  fraction or mantissa \code{f} and integer exponent \code{E}, it returns
  \eqn{x = f * 2^E}.
  Viewed differently, it's the fastest way to multiply or divide (double
  precision) numbers with \eqn{2^E}.
}
\usage{
frexp(x)
ldexp(f, E)
}
\arguments{
  \item{x}{numeric (coerced to \code{double}) vector.}
  \item{f}{numeric fraction (vector), in \eqn{[0.5, 1)}.}
  \item{E}{integer valued, exponent of \code{2}, i.e., typically in
    \code{(-1024-50):1024}, otherwise the result will underflow to 0 or
    overflow to \code{+/- Inf}.}
}
%% \details{
%% }
\value{
  \code{frexp} returns a \code{\link{list}} with named components \code{r}
  (of type \code{double}) and \code{e} (of type \code{integer}).
}
\references{
  On unix-alikes, typically \command{man frexp} and \command{man ldexp}
}
\author{Martin Maechler}
%% \note{
%% }
\seealso{
  Vaguely relatedly, \code{\link{log1mexp}()}, \code{\link{lsum}}, \code{\link{logspace.add}}.
}
\examples{
set.seed(47)
x <- c(0, 2^(-3:3), (-1:1)/0,
       rlnorm(2^12, 10, 20) * sample(c(-1,1), 512, replace=TRUE))
head(x, 11)
rF <- frexp(x)
sapply(rF, summary) # (nice only when x had no NA's ..)
ar <- abs(rF$r)
stopifnot(0.5 <= ar[is.finite(x) & x != 0], ar[is.finite(x)] < 1,
          is.integer(rF$e))
ldx <- with(rF, ldexp(r, e))
stopifnot(all.equal(x, ldx, tol= 4*.Machine$double.eps))
## but actually, they should even be identical :
identical(x, ldx) # TRUE  on all (? / decent) platforms
}
\keyword{arith}

