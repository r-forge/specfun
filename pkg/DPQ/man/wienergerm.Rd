\name{pchisqW.R}
\alias{pchisqW.R}
\title{Wienergerm Approximations to pchisq() the Chisquare Distribution}
\description{
  Compute Wienergerm approximations to \code{\link{pchisq}()}, the
  (non-central) chisquare distribution.
}
\usage{
pchisqW.R(x, df, ncp = 0, lower.tail = TRUE, log.p = FALSE,
          variant = c("s", "f"), verbose = getOption("verbose"))
}

\arguments{
  \item{x,df,ncp}{main arguments, see  \code{\link{pchisq}}.}
  \item{lower.tail,log.p}{\code{\link{logical}}, see \code{\link{pchisq}}.}

  \item{variant}{a string, currently either \code{"f"} for the \bold{f}irst or
    \code{"s"} for the \bold{s}econd Wienergerm approximation in
    Penev \& Raykov (2000).}
  \item{verbose}{logical/integer ...}% FIXME
}
\value{
  bla FIXME
}
\author{Martin Maechler}
%% \details{   TODO -- FIXME
%%
%% ==> info in header of ../src/wienergerm_nchisq.c
%% but also in Martin's extensive example scripts, and hand notes
%% }
\references{
  Penev, Spiridon and Raykov, Tenko (2000)
  A Wiener Germ approximation of the noncentral chi square
  distribution and of its quantiles.
  Computational Statistics \bold{15}, 219--228.
  \doi{10.1007/s001800000029}% ~/save/papers/Numerics/a-wiener-germ-approximation.ps

  Dinges, H. (1989)
  Special cases of second order Wiener germ approximations.
  \emph{Probability Theory and Related Fields}, \bold{83}, 5--57.
}
\seealso{
  \code{\link{pchisq}}
}
\examples{
"TODO"
}
\keyword{distribution}
\keyword{math}
