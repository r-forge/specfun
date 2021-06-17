\name{dgamma-utils}
\title{Utility Functions for \code{dgamma()} -- Pure R Versions}
\alias{bd0}
\alias{bd0_p1l1d1}
\alias{bd0_p1l1d}
\alias{bd0_l1pm}
\alias{bd0C}
\alias{ebd0}
\alias{ebd0C}
\alias{dpois_raw}
\alias{stirlerr}
\alias{lgammacor}
\description{
  Mostly, pure \R transcriptions of the C code utility functions for
  \code{\link{dgamma}()} and similar \dQuote{base} density functions by
  Catherine Loader.

  \code{bd0C()} interfaces to C code which corresponds to \R's C Mathlib
  \code{bd0()}.

  These have extra arguments with defaults that correspond
  to the C Mathlib code hardwired cutoffs and tolerances.
}
\usage{% >>> ../R/dgamma.R <<<
dpois_raw(x, lambda, log=FALSE,
          version, % =  c("bd0_v1", "bd0_p1l1d", "bd0_p1l1d1", "bd0_l1pm", "ebd0_v1")
          ## the defaults for version will probably change in the future
          bd0.delta = 0.1,
          ## optional arguments of log1pmx() :
          tol_logcf = 1e-14, eps2 = 0.01, minL1 = -0.79149064, trace.lcf = verbose,
          logCF = if (is.numeric(x)) logcf else logcfR,
          verbose = FALSE)
bd0(x, np,
    delta = 0.1, maxit = as.integer(-1100 / log2(delta)),
    s0 = .Machine$double.xmin,
    verbose = getOption("verbose"))
bd0C(x, np, delta = 0.1, maxit = 1000L, version = "R4.0", verbose = getOption("verbose"))
# "simple" log1pmx() based versions :
bd0_p1l1d1(x, M, tol_logcf = 1e-14, ...)
bd0_p1l1d (x, M, tol_logcf = 1e-14, ...)
bd0_l1pm  (x, M, tol_logcf = 1e-14, ...)

ebd0 (x, M, verbose = getOption("verbose"), ...) # experimental, may disappear !!
ebd0C(x, M, verbose = getOption("verbose"))

stirlerr(n, scheme = c("R3", "R4.1"),
         cutoffs = switch(scheme
                        , R3   = c(15, 35, 80, 500)
                        , R4.1 = c(7.5, 8.5, 10.625, 12.125, 20, 26, 55, 200, 3300)
                          ),
         use.halves = missing(cutoffs),
         verbose = FALSE)
lgammacor(x, nalgm = 5, xbig = 2^26.5)% --> ../R/utils.R <-> ../src/lgammacor.c
}
\arguments{
  \item{x, n}{\code{\link{numeric}} (or number-alike such as "mpfr").}
  \item{lambda, np, M}{each \code{\link{numeric}} (or number-alike ..); distribution parameters.}
  \item{log}{logical indicating if the log-density should be returned,
    otherwise the density at \code{x}.}
  \item{verbose}{logical indicating if some information about the
    computations are to be printed.}
  \item{delta, bd0.delta}{a positive number \eqn{< 1} (practically required
    to be \eqn{\le .99}), a cutoff for \code{bd0()}
    where the \code{\link{logcf}()} series expansion is used when
    \eqn{|x - M| < delta*(x+M)}.}
  \item{tol_logcf, eps2, minL1, trace.lcf, logCF, ...}{
    optional tuning arguments passed to \code{\link{log1pmx}()}.}
  \item{maxit}{the number of \code{logcf()} terms to be used in
    \code{bd0()} when \eqn{|x-M|} is small.  The default is \eqn{k} such
    that \eqn{\delta^{2k} \le 2^{-1022-52}}, i.e., will underflow to zero.}
  \item{s0}{the very small \eqn{s_0} determining that \code{bd0() = s}
    already before the locf series expansion.}
  \item{version}{a \code{\link{character}} string specifying the version of
    \eqn{bd0()} used.}
  \item{scheme}{a \code{\link{character}} string specifying the
    \code{cutoffs} scheme.}
  \item{cutoffs}{an increasing numeric vector, required to start with
    with \code{cutoffs[1] <= 15} specifying the cutoffs to switch from 2 to 3 to ...,
    up to 10 term approximations for non-small \code{n}, where the direct
    formula loses precision.  When missing (as by default), \code{scheme}
    is used, where \code{scheme = "R3"} chooses (15, 35, 80, 500), the
    cutoffs in use in \R versions up to (and including) 4.0.z.}
  \item{use.halves}{\code{\link{logical}} indicating if the full-accuracy
    prestored values should be use when
    \eqn{2n \in \{0,1,\dots,30\}}{2n in {0,1,..,30}}, i.e.,
    \eqn{n \le 15}{n <= 15} and n is integer or integer +
    \eqn{\frac{1}{2}}{1/2}.
    Turn this off to judge the underlying approximation accuracy by
    comparison with MPFR.   However, keep the default \code{TRUE} for
    back-compatibility.}
  \item{nalgm}{number of terms to use for Chebyshev polynomial approxmation
    in \code{lgammacor()}.  The default, 5, is the value hard wired in \R's
    C Mathlib.}
  \item{xbig}{a large positive number; if \code{x >= xbig}, the simple
    asymptotic approximation \code{lgammacor(x) := 1/(12*x)} is used.  The
    default, \eqn{2^{26.5} = 94906265.6}, is the value hard wired in \R's C
    Mathlib.}
}
\value{
  a numeric vector \dQuote{like} \code{x}; in some cases may also be an
  (high accuracy) "mpfr"-number vector, using CRAN package \CRANpkg{Rmpfr}.

  \code{lgammacor(x)} originally returned \code{NaN} for all \eqn{|x| < 10},
  as its Chebyshev polynomial approximation has been constructed for
  \eqn{x \in [10, xbig]}{x in [10, xbig]},
  specifically for \eqn{u \in [-1,1]} where
  \eqn{t := 10/x \in [1/x_B, 1]}  and
  \eqn{u := 2t^2 -1 \in [-1 + \epsilon_B, 1]}.
}
\details{
  \describe{% >> TODO: Write Vignette w/ proper LaTeX :
	    % ===>  ../vignettes/log1pmx-etc.Rnw
    \item{\code{bd0()}:}{Loader's \dQuote{Binomial Deviance} function; for
      \eqn{x, M > 0} (where the limit \eqn{x \to 0}{x -> 0} is allowed).
      In the case of \code{\link{dbinom}}, \eqn{x} are integers (and
      \eqn{M = n p}), but in general \code{x} is real.

      \deqn{bd_0(x,M) := M \cdot D_0\bigl(\frac{x}{M}\bigr),}{bd_0(x,M) := M * D_0(x/M),} where
      \eqn{D_0(u) := u \log(u) + 1-u = u(\log(u) - 1) + 1}.  Hence
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
      the Taylor series expansion is useful for small \eqn{|t|}.
    }% item
  }% describe
}
\references{
  C. Loader (2000), see \code{\link{dbinom}}'s documentation.
  % ~/R/D/r-devel/R/src/library/stats/man/Binomial.Rd
}
\author{Martin Maechler}
\seealso{
  \code{\link{dgamma}}, % ~/R/D/r-devel/R/src/library/stats/man/GammaDist.Rd
  \code{\link{dpois}}.
  High precision versions \code{\link[DPQmpfr]{stirlerrM}(n)} and
  \code{stirlerrSer(n,k)} in package \CRANpkg{DPQmpfr} (via the
  \CRANpkg{Rmpfr} and \CRANpkg{gmp} packages).
}
\examples{
n <- seq(1, 50, by=1/4)
st.n <- stirlerr(n) # now vectorized
stopifnot(identical(st.n, sapply(n, stirlerr)))
plot(n, st.n, type = "b", log="xy", ylab = "stirlerr(n)")

x <- 800:1200
bd0x1k <- bd0(x, np = 1000)
plot(x, bd0x1k, type="l", ylab = "bd0(x, np=1000)")
bd0x1kC <- bd0C(x, np = 1000)
lines(x, bd0x1kC, col=2)
bd0.1d1 <- bd0_p1l1d1(x, 1000)
bd0.1d  <- bd0_p1l1d (x, 1000)
bd0.1pm <- bd0_l1pm  (x, 1000)
stopifnot(exprs = {
    all.equal(bd0x1kC, bd0x1k,  tol=1e-15) # even tol=0 currently ..
    all.equal(bd0x1kC, bd0.1d1, tol=1e-15)
    all.equal(bd0x1kC, bd0.1d , tol=1e-15)
    all.equal(bd0x1kC, bd0.1pm, tol=1e-15)
})

str(log1pmx) ##--> play with  { tol_logcf, eps2, minL1, trace.lcf, logCF }

ebd0x1k <- ebd0 (x, 1000) # FIXME: subscript out of bout
exC     <- ebd0C(x, 1000)
stopifnot(all.equal(exC, ebd0x1k, tol=4e-16))
lines(x, colSums(ebd0x1k), col=adjustcolor(4, 1/2), lwd=4)

x <- 0:250
dp   <- dpois    (x, 48, log=TRUE)# R's 'stats' pkg function
dp.r <- dpois_raw(x, 48, log=TRUE)
all.equal(dp, dp.r, tol = 0) # on Linux 64b, see TRUE
stopifnot(all.equal(dp, dp.r, tol = 1e-14))
## dpois_raw()  versions:
(vers <- eval(formals(dpois_raw)$version))
mv <- sapply(vers, function(v) dpois_raw(x, 48, version=v))
matplot(x, mv, type="h", log="y") # "fine"
if(all(abs(1 - mv[,"ebd0_C1"]/mv[,"ebd0_v1"]) < 2e-16)) mv <- mv[, vers != "ebd0_C1"]
## now look at *relative* errors -- need "Rmpfr" for "truth"
if(requireNamespace("Rmpfr")) {

    dM <- Rmpfr::dpois(Rmpfr::mpfr(x, 256), 48)
    asN <- Rmpfr::asNumeric
    relE <- asN(mv / dM - 1)
    cols <- adjustcolor(1:ncol(mv), 1/2)

    mtit <- "relative Errors of dpois_raw(x, 48, version = * )"
    matplot(x, relE, type="l", col=cols, lwd=3, lty=1, main=mtit)
    legend("topleft", colnames(mv), col=cols, lwd=3, bty="n")

    matplot(x, abs(relE), ylim=pmax(1e-18, range(abs(relE))), type="l", log="y",
            main=mtit, col=cols, lwd=2, lty=1, yaxt="n")
    sfsmisc::eaxis(2)
    legend("bottomright", colnames(mv), col=cols, lwd=2, bty="n", ncol=3)
    ee <- c(.5, 1, 2)* 2^-52; eC <- quote(epsilon[C])
    abline(h=ee, lty=2, col="gray", lwd=c(1,2,1))
    axis(4, at=ee[2:3], expression(epsilon[C], 2 * epsilon[C]), col="gray", las=1)
    par(new=TRUE)
    plot(x, asN(dM), type="h", col=adjustcolor("darkgreen", 1/3), axes=FALSE, ann=FALSE)
    stopifnot(abs(relE) < 8e-13) # seen 2.57e-13
}# Rmpfr
}
\keyword{distribution}
