\name{dpsifn}
\alias{dpsifn}
\title{Psi Gamma Functions Workhorse from R's API}
\description{
  Log Gamma derivatives,  Psi Gamma functions.  \code{dpsifn()} is an \R
  interface to the \R API function \code{R_dpsifn()}.
}
\usage{
dpsifn(x, m, deriv1 = 0L, k2 = FALSE)
}
\arguments{
  \item{x}{numeric vector.}
  \item{m}{number of derivatives to return, an integer >= 0.}
  \item{deriv1}{\dQuote{start} derivative ....}%FIXME
  \item{k2}{a \code{\link{logical}} specifying if \code{kode = 2} should be
    applied.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A \code{\link{list}} with three named components,
  \item{dpsi}{numeric matrix of dpshi() values.}
  \item{overflow}{Description of 'comp1'}% FIXME
  \item{ierr}{Description of 'comp2'}% FIXME
}
\references{
  See those in \code{\link{psigamma}}
}
\author{Martin Maechler (R interface); R Core et al., see
  \code{\link{digamma}}.
}
%% \note{
%% }
\seealso{
  \code{\link{digamma}}, \code{\link{trigamma}}, \code{\link{psigamma}}.
}
\examples{
x <- seq(-5,3, by=1/4)
dpx <- dpsifn(x, m=5)# default deriv = 0
dpx # in R <= 4.2.1, see that sometimes the 'nz' (under-over-flow count) was uninitialized !!
}
\keyword{math}
\concept{psi gamma functions}
