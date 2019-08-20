\name{qchisqAppr}
\title{Compute Approximate Quantiles of the Chi-Squared Distribution}
\alias{qchisqKG}
\alias{qchisqWH}
\alias{qchisqAppr}
\alias{qchisqAppr.R}
\description{
  Compute quantiles (inverse distribution values) for the chi-squared distribution.
  using Johnson,Kotz,.. ............TODO.......
}
\usage{
%% *central* chi^2 : R code in ../R/qgamma-fn.R  : qchisqAppr.R()  & ../R/qnchisq.R
qchisqKG    (p, df, lower.tail = TRUE, log.p = FALSE)
qchisqWH    (p, df, lower.tail = TRUE, log.p = FALSE)
qchisqAppr  (p, df, lower.tail = TRUE, log.p = FALSE, tol = 5e-7)
qchisqAppr.R(p, df, lower.tail = TRUE, log.p = FALSE, tol = 5e-07,
             maxit = 1000, verbose = getOption("verbose"), kind = NULL)
}
\arguments{%--> ~/R/D/r-devel/R/src/library/stats/man/Chisquare.Rd :
%% ===>  qchisq(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  \item{p}{vector of probabilities.}
  \item{df}{degrees of freedom \eqn{> 0}, maybe non-integer; must have
    length 1.}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qchisq}()}; must
  have length 1.}
%%% TODO: \item{method}{a string specifying the approximation method to be used.}
%%% -----  INSTEAD of all these functions ??
  \item{tol}{non-negative number, the convergence tolerance}
  %% \item{g}{??? FIXME eliminate ?????? FIXME ?????}
  \item{maxit}{the maximal number of iterations}
  \item{verbose}{logical indicating if the algorithm should produce
    \dQuote{monitoring} information.}
  \item{kind}{the \emph{kind} of approximation; if \code{NULL}, the
    default, the approximation chosen depends on the arguments; notably it
    is chosen separately for each \code{p}.  Otherwise, it must be a
    \code{\link{character}} string.  The main approximations are
    Wilson-Hilferty versions, when the string contains \code{"WH"}.
    More specifically, it must be one of the strings
    \describe{
      \item{"chi.small"}{particularly useful for small chisquared values \code{p};... ...}
      \item{"WH"}{... ...}
      \item{"p1WH"}{... ...}
      \item{"WHchk"}{... ...}
      \item{"df.small"}{particularly useful for small degrees of freedom \code{df}... ...}
    }
  }
}
%% \details{
%% }
\value{
  ...
}
%% \references{
%%
%% }
\author{Martin Maechler}
\seealso{
  \code{\link{qchisq}}.   Further, our approximations to the
  \emph{non-central} chi-squared quantiles, \code{\link{qnchisqAppr}}
}
\examples{
  ## TODO
}
\keyword{distribution}
