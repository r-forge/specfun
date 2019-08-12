\name{lbeta}
\alias{lbetaM}
\alias{lbeta.asy}
\alias{lbetaMM}
\alias{lbeta.n}
\alias{lbeta.n.half}
\title{Log Beta Approximations}
\description{
  Compute \code{log(\link{beta}(a,b))} in a simple (fast) or asymptotic way.
}
\usage{%% in ../R/beta-fns.R <<<<<<<
lbetaM   (a, b, k.max = 5, give.all = FALSE)
lbeta.asy(a, b, k.max = 5, give.all = FALSE)
lbetaMM  (a, b, Cut.asy = 1e-2, verbose = FALSE)

lbeta.n     (n, a)
lbeta.n.half(n, a)
}
\arguments{
  \item{a, b}{the Beta parameters, see \code{\link{beta}}.}
  \item{k.max}{ .. }
  \item{give.all}{\code{\link{logical}} .. }
  \item{Cut.asy}{cutoff value from where to switch to asymptotic formula.}%% FIXME details
  \item{verbose}{logical (or integer) indicating if and how much monitoring
    information should be printed to the console.}
  \item{n}{ .. }
}
%% \details{
%% }
\value{
  ..
}
%% \references{
%% }
\author{Martin Maechler}
%% \note{
%% }
\seealso{
  \R's \code{\link{beta}} function
}
\examples{
  ## TODO
}
\keyword{math}

