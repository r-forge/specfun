\name{dgamma-utils}
\title{Binomial Deviance -- Auxiliary Functions for \code{dgamma()} Etc}
\alias{bd0}
\alias{bd0_p1l1d1}
\alias{bd0_p1l1d}
\alias{bd0_p1l1}
\alias{bd0_l1pm}
\alias{bd0C}
\alias{ebd0}
\alias{ebd0C}
\alias{dpois_raw}
\alias{dpois_simpl}
\alias{dpois_simpl0}
\description{
  The \dQuote{binomial deviance} function
  \code{bd0(x, M) :=} \eqn{D_0(x,M) := M \cdot d_0(x/M)}{D0(x,M) := M * d_0(x/M)},
  where \eqn{d_0(r) := r\log(r) + 1-r \stackrel{!}{=} p_1l_1(r-1)}{d_0(r) := r*log(r) + 1-r = p1l1(r-1)}.
  Mostly, pure \R transcriptions of the C code utility functions for
  \code{\link{dgamma}()}, \code{\link{dbinom}()}, \code{\link{dpois}()}, \code{\link{dt}()},
  and similar \dQuote{base} density functions by Catherine Loader. \cr
  These have extra arguments with defaults that correspond
  to \R's Mathlib C code hardwired cutoffs and tolerances.
}% .. more in \details below

\usage{% >>> ../R/dgamma.R <<<
dpois_raw(x, lambda, log=FALSE,
          version, % =  c("bd0_v1", "bd0_p1l1d", "bd0_p1l1d1", "bd0_p1l1", "bd0_l1pm", "ebd0_v1")
          ## the defaults for `version` will probably change in the future
          small.x__lambda = .Machine$double.eps,
          bd0.delta = 0.1,
          ## optional arguments of log1pmx() :
          tol_logcf = 1e-14, eps2 = 0.01, minL1 = -0.79149064, trace.lcf = verbose,
          logCF = if (is.numeric(x)) logcf else logcfR,
          verbose = FALSE)

dpois_simpl (x, lambda, log=FALSE)
dpois_simpl0(x, lambda, log=FALSE)

bd0(x, np, delta = 0.1,
    maxit = max(1L, as.integer(-1100 / log2(delta))),
    s0 = .Machine$double.xmin,
    verbose = getOption("verbose"))
bd0C(x, np, delta = 0.1, maxit = 1000L, version = c("R4.0", "R_2025_0510"),
     verbose = getOption("verbose"))
# "simple" log1pmx() based versions :
bd0_p1l1d1(x, M, tol_logcf = 1e-14, ...)
bd0_p1l1d (x, M, tol_logcf = 1e-14, ...)
# p1l1() based version {possibly using log1pmx(), too}:
bd0_p1l1  (x, M, ...)
# using faster formula for large |M-x|/x
bd0_l1pm  (x, M, tol_logcf = 1e-14, ...)

ebd0 (x, M, verbose = getOption("verbose"), ...) # experimental, may disappear !!
ebd0C(x, M, verbose = getOption("verbose"))
}
\arguments{
  \item{x}{\code{\link{numeric}} (or number-alike such as "mpfr").}
  \item{lambda, np, M}{each \code{\link{numeric}} (or number-alike ..); distribution parameters.}
  \item{log}{logical indicating if the log-density should be returned,
    otherwise the density at \code{x}.}
  \item{verbose}{logical indicating if some information about the
    computations are to be printed.}
  \item{small.x__lambda}{small positive number; for \code{dpois_raw(x, lambda)},
    when \code{x/lambda} is not larger than \code{small.x__lambda}, the
    direct log poisson formula is used instead of \code{ebd0()} or
    \code{bd0()} and \code{\link{stirlerr}()}.}
  \item{delta, bd0.delta}{a non-negative number \eqn{\delta < 1}, a cutoff
    for \code{bd0()} where a continued fraction series expansion is used
    when \eqn{|x - M| \le \delta\cdot(x+M)}{|x - M| <= delta*(x+M)}.}
  \item{tol_logcf, eps2, minL1, trace.lcf, logCF, ...}{
    optional tuning arguments passed to \code{\link{log1pmx}()}, and to its
    options passed to \code{\link{logcf}()}.}
  \item{maxit}{the number of series expansion terms to be used in
    \code{bd0()} when \eqn{|x-M|} is small.  The default is \eqn{k} such
    that \eqn{\delta^{2k} \le 2^{-1022-52}}, i.e., will underflow to zero.}
  \item{s0}{the very small \eqn{s_0} determining that \code{bd0() = s}
    already before the locf series expansion.}
  \item{version}{a \code{\link{character}} string specifying the version of
    \eqn{bd0()} to use.  The versions available \emph{and} the current
    default are partly experimental!}
}
\details{
  \describe{% >> Now, most of this and more is in vignette
	    % ../vignettes/log1pmx-etc.Rnw   <<<<<<<<<
    \item{\code{bd0()}:}{Loader's \dQuote{Binomial Deviance} function; for
      \eqn{x, M > 0} (where the limit \eqn{x \to 0}{x -> 0} is allowed).
      In the case of \code{\link{dbinom}}, \eqn{x} are integers (and
      \eqn{M = n p}), but in general \code{x} is real.

      \deqn{bd_0(x,M) := M \cdot D_0\bigl(\frac{x}{M}\bigr),}{bd_0(x,M) := M * D_0(x/M),} where
      \eqn{D_0(u) := u \log(u) + 1-u = u(\log(u) - 1) + 1.}{D_0(u) := u log(u) + 1-u = u(log(u) - 1) + 1.}  Hence
      \deqn{bd_0(x,M) = M \cdot \bigl(\frac{x}{M}(\log(\frac{x}{M}) -1) +1 \bigr) =
	x \log(\frac{x}{M}) - x + M.}{bd_0(x,M) = M *((x/M)*(log(x/M) - 1) +1) = x log(x/M) - x + M.}

      A different way to rewrite this from Martyn Plummer, notably for important situation when
      \eqn{\left|x-M \right| \ll M}{|x-M| << M}, is using \eqn{t := (x-M)/M}
      (and \eqn{\left|t \right| \ll 1}{|t| << 1} for that situation),
      equivalently, \eqn{\frac{x}{M} = 1+t}{x/M = 1+t}.
      Using \eqn{t},
      \deqn{bd_0(x,M) = \log(1+t) - t \cdot M = M \cdot [(t+1)(\log(1+t) - 1) + 1]
                                              = M \cdot [(t+1) \log(1+t) - t]
                                              = M \cdot p_1l_1(t),}{
            bd_0(x,M) = log(1+t) - t * M = M * ((t+1)(log(1+t) - 1) + 1)
                                         = M * ((t+1) log(1+t) - t) = M * p1l1(t),}
     and \deqn{p_1l_1(t) := (t+1)\log(1+t) - t = \frac{t^2}{2} - \frac{t^3}{6} ...}{
                p1l1 (t) := (t+1)*log(1+t) - t = t^2/2 - t^3/6 ...} where
      the Taylor series expansion is useful for small \eqn{|t|}, see \code{\link{p1l1}}.

      Note that \code{bd0(x, M)} now also works when \code{x} and/or
      \code{M} are arbitrary-accurate mpfr-numbers (package \CRANpkg{Rmpfr}).

      \code{bd0C()} interfaces to C code which corresponds to \R's C Mathlib (Rmath) \code{bd0()}.
    }% bd0
  }% describe
}
\value{
  a numeric vector \dQuote{like} \code{x}; in some cases may also be an
  (high accuracy) "mpfr"-number vector, using CRAN package \CRANpkg{Rmpfr}.

  \code{ebd0()} (\R code) and \code{ebd0C()} (interface to \command{C}
  code) are \emph{experimental}, meant to be precision-extended version of
  \code{bd0()}, returning \code{(yh, yl)} (high- and low-part of \code{y},
  the numeric result).  In order to work for \emph{long} vectors \code{x},
  \code{yh, yl} need to be \code{\link{list}} components; hence we return a
  two-column \code{\link{data.frame}} with column names \code{"yh"} and
  \code{"yl"}.
}
\references{
  C. Loader (2000), see \code{\link{dbinom}}'s documentation.
  % ~/R/D/r-devel/R/src/library/stats/man/Binomial.Rd

  Martin \enc{Mächler}{Maechler} (2021 ff)
  %log1pmx(), bd0(), stirlerr() -- Computing Poisson, Binomial, Gamma Probabilities in \R}
  \emph{\code{log1pmx}, \dots Computing \dots Probabilities in \R}.
  \CRANpkg{DPQ} package vignette % ../vignettes/log1pmx-etc.Rnw
  \url{https://CRAN.R-project.org/package=DPQ/vignettes/log1pmx-etc.pdf}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{p1l1}()} and its variants, e.g., \code{p1l1ser()};
  \code{\link{log1pmx}()} which is called from \code{bd0_p1l1d*()} and
  \code{bd0_l1pm()}.  Then,
  \code{\link{stirlerr}} for Stirling's error function, % -- ./stirlerr.Rd
  complementing \code{bd0()} for computation of Gamma, Beta, Binomial and Poisson probabilities.
  \R's own
  \code{\link{dgamma}}, % ~/R/D/r-devel/R/src/library/stats/man/GammaDist.Rd
  \code{\link{dpois}}.
}
\examples{
x <- 800:1200
bd0x1k <- bd0(x, np = 1000)
plot(x, bd0x1k, type="l", ylab = "bd0(x, np=1000)")
bd0x1kC <- bd0C(x, np = 1000)
lines(x, bd0x1kC, col=2)
bd0.1d1 <- bd0_p1l1d1(x, 1000)
bd0.1d  <- bd0_p1l1d (x, 1000)
bd0.1p1 <- bd0_p1l1  (x, 1000)
bd0.1pm <- bd0_l1pm  (x, 1000)
stopifnot(exprs = {
    all.equal(bd0x1kC, bd0x1k,  tol=1e-14) # even tol=0 currently ..
    all.equal(bd0x1kC, bd0.1d1, tol=1e-14)
    all.equal(bd0x1kC, bd0.1d , tol=1e-14)
    all.equal(bd0x1kC, bd0.1p1, tol=1e-14)
    all.equal(bd0x1kC, bd0.1pm, tol=1e-14)
})

str(log1pmx) #--> play with { tol_logcf, eps2, minL1, trace.lcf, logCF }  --> vignette

ebd0x1k <- ebd0 (x, 1000)
exC     <- ebd0C(x, 1000)
stopifnot(all.equal(exC, ebd0x1k,          tol = 4e-16),
          all.equal(bd0x1kC, rowSums(exC), tol = 2e-15))
lines(x, rowSums(ebd0x1k), col=adjustcolor(4, 1/2), lwd=4)

## very large args --> need new version "R_2025..."
np <- 117e306;  np / .Machine$double.xmax # 0.65
x <- (1:116)*1e306
tail(bd0L <-  bd0C(x, np, version = "R4")) # underflow to 0
tail(bd0Ln <- bd0C(x, np, version = "R_2025_0510"))
bd0L.R <- bd0(x, np)
all.equal(bd0Ln, bd0L.R, tolerance = 0) # see TRUE
stopifnot(exprs = {
    is.finite(bd0Ln)
    bd0Ln > 4e303
    tail(bd0L, 54) == 0
    all.equal(bd0Ln, np*p1l1((x-np)/np), tolerance = 5e-16)
    all.equal(bd0Ln, bd0L.R, tolerance = 5e-16)
})

%% --> ../vignettes/log1pmx-etc.Rnw
}
\keyword{arith}
\keyword{distribution}
