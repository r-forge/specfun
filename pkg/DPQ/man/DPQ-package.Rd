\name{DPQ-package}
\alias{DPQ-package}
\alias{DPQ}
\docType{package}
\title{\packageTitle{DPQ}}
\description{
  \packageDescription{DPQ}
}
\details{
  The DESCRIPTION file: % >> ../DESCRIPTION <<
  \packageDESCRIPTION{DPQ}
  \packageIndices{DPQ}

  An important goal is to investigate diverse algorithms and approximations
  of \R's own density (\code{d*()}), probability (\code{p*()}), and
  quantile (\code{q*()}) functions, notably in \dQuote{border} cases where
  the traditional published algorithms have shown to be suboptimal, not
  quite accurate, or even useless.

  Examples are border cases of the beta distribution, or \bold{non-central}
  distributions such as the non-central chi-squared and t-distributions.
}
\author{
  Principal author and maintainer: \packageMaintainer{DPQ}
  %% \packageAuthor{DPQ}
}
%% \references{
%% ~~ Literature or other references for background information ~~
%% }
\seealso{
  The package \pkg{\link[DPQmpfr:DPQmpfr-package]{DPQmpfr}},
  which builds on this package and on
  \pkg{\link[Rmpfr:Rmpfr-package]{Rmpfr}}.
}
\examples{
\donttest{
## Show problem in R's non-central t-distrib. density (and approximations):
example(dntJKBf)
}
%% TODO add other few cases (be *fast*!), maybe  non-central chisq
}
\keyword{ distribution }
\keyword{ math }
\keyword{ package }
