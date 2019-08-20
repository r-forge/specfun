\name{pnchisqWienergerm}
\title{Wienergerm Approximations to (Non-Central) Chisquare Probabilities}
%% TODO: also read/use >>>>  ../wienergerm_README <<<<<<<<<<<<<<<
\alias{pchisqV}% w/o "n" -- NB:  --> ./pnchisqAppr.Rd for 'pnchisqV() % *with* "n"
\alias{pchisqW}
\alias{pchisqW.}
\alias{pchisqW.R}
\alias{sW}
\alias{qs}
\alias{z0}
\alias{z.f}
\alias{z.s}
%-
\alias{g2}
\alias{gnt}
\alias{h}
\alias{h0}
\alias{h1}
\alias{h2}
\alias{hnt}
\alias{scalefactor}% ? no longer exported
\alias{ss}
\alias{ss2}
\alias{ss2.}
\description{%%--> ../R/wienergerm_nchisq-fn.R <--

  Functions implementing the two Wiener germ approximations to
  \code{\link{pchisq}()}, the (non-central) chisquare distribution, and to
  \code{\link{qchisq}()} its inverse, the quantile function.

  These have been proposed by Penev and Raykov (2000) who also listed a
  Fortran implementation.

  In order to use them in numeric boundary cases, Martin Maechler has
  improved the original formulas.

  \bold{Auxiliary functions:}
  \describe{
    \item{\code{sW()}:}{The \eqn{s()} as in the Wienergerm approximation,
    but using Taylor expansion when needed, i.e., \code{(x*ncp / df^2) << 1}.}

    \item{\code{qs()}:}{ ... }
    \item{\code{z0()}:}{ ... }
    \item{\code{z.f()}:}{ ... }
    \item{\code{z.s()}:}{ ... }

    ..................
    ..................

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

%%------

ss   (x, df, ncp, i.max = 10000)
ss2  (x, df, ncp, i.max = 10000, eps = .Machine$double.eps)
ss2. (q, df, ncp = 0, errmax = 1e-12, reltol = 2 * .Machine$double.eps,
      itrmax = 1e+05, eps = reltol)

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
  \item{verbose}{logical (or integer) indicating if or how much diagnostic
    output should be printed to the console during the computations.}
  \item{f.s}{ .... }
  \item{eps1}{ .... }
  \item{sMax}{ .... }

  \item{i.max}{ .. }
  \item{eps}{ .. }
  \item{errmax}{ .. }
  \item{reltol}{ .. }
  \item{itrmax}{ .. }
}
\note{
  The exact auxiliary function names etc, are still considered
  \emph{provisional}; currently they are exported for easier documentation
  and use, but may well all disappear from the exported functions or even
  completely.
}
\details{
  ....TODO... or write vignette
%% ==> info in header of ../src/wienergerm_nchisq.c
%% but also in Martin's extensive example scripts, and hand notes
}
\value{
  bla FIXME
}
\author{Martin Maechler, mostly end of Jan 2004}% <--> ../wienergerm_README <<<<
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
% ../tests/wienergerm-accuracy.R
% ../tests/wienergerm-pchisq-tst.R
% ../tests/wienergerm_nchisq.R
}
\keyword{distribution}
\keyword{math}
