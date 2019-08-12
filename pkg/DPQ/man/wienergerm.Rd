\name{pchisqW}
\title{Wienergerm Approximations to pchisq() the Chisquare Distribution}
\alias{pchisqV}
\alias{pchisqW}
\alias{pchisqW.}
\alias{pchisqW.R}
\alias{sW}
\alias{qs}
\alias{z0}
\alias{z.f}
\alias{z.s}
\description{
  Compute Wienergerm approximations to \code{\link{pchisq}()}, the
  (non-central) chisquare distribution.

  \bold{Auxiliary functions:}
  \describe{
    \item{\code{sW()}:}{The \eqn{s()} as in the Wienergerm approximation,
    but using Taylor expansion when needed, i.e., \code{(x*ncp / df^2) << 1}.}

    \item{\code{qs()}:}{ ... }
    \item{\code{z0()}:}{ ... }
    \item{\code{z.f()}:}{ ... }
    \item{\code{z.s()}:}{ ... }
  } % descr.
}
\usage{
% R code in  ../R/wienergerm_nchisq-fn.R :
%                 ----------------------
pchisqW. (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE,
          Fortran = TRUE, variant = c("s", "f"))
pchisqV  (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE,
          Fortran = TRUE, variant = c("s", "f"))
pchisqW  (q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE, variant = c("s", "f"))
pchisqW.R(x, df, ncp = 0, lower.tail = TRUE, log.p = FALSE, variant = c("s", "f"),
          verbose = getOption("verbose"))

sW(x, df, ncp)
qs(x, df, ncp, f.s = sW(x, df, ncp), eps1 = 1/2, sMax = 1e+100)
z0(x, df, ncp)
z.f(x, df, ncp)
z.s(x, df, ncp, verbose = getOption("verbose"))
}

\arguments{
  \item{q,x}{vector of quantiles (main argument, see  \code{\link{pchisq}}).}
  \item{df}{degrees of freedom (non-negative, but can be non-integer).}
  \item{ncp}{non-centrality parameter (non-negative).}
  \item{lower.tail,log.p}{\code{\link{logical}}, see \code{\link{pchisq}}.}

  \item{variant}{a \code{\link{character}} string, currently either
    \code{"f"} for the \bold{f}irst or
    \code{"s"} for the \bold{s}econd Wienergerm approximation in
    Penev \& Raykov (2000).}
  \item{Fortran}{logical specifying if the Fortran or the C version should
    be used.}
  \item{verbose}{logical/integer ...}% FIXME
  \item{f.s}{ .... }
  \item{eps1}{ .... }
  \item{sMax}{ .... }
}
\value{
  bla FIXME
}
\author{Martin Maechler, 26 Jan 2004}

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
  \code{\link{pchisq}}, and other approximations for it:
  \code{\link{pnchisq}()} etc.
}
\examples{
## see  example(pnchisqAppr)   which looks at all of the pchisq() approximating functions
}
\keyword{distribution}
\keyword{math}
