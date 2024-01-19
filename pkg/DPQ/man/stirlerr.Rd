\name{stirlerr}%%  was part of >>> ./dgamma-utils.Rd <<<
\title{Stirling's Error Function - Auxiliary for Gamma, Beta, etc}% Binomial, Poisson Distributions
\alias{stirlerr}
\alias{stirlerrC}
\alias{stirlerr_simpl}
\alias{lgammacor}
\description{%% -- ~ = ~ as ../../DPQmpfr/man/stirlerrM.Rd
  Stirling's approximation (to the factorial or \eqn{\Gamma} function)
  error in \eqn{\log}{log} scale is the difference of the left and right hand
  side of Stirling's approximation to \eqn{n!},%
  \eqn{n! \approx \bigl(\frac{n}{e}\bigr)^n \sqrt{2\pi n},}{%
       n! ~= (n/e)^n * sqrt(2*pi*n),} i.e., \code{stirlerr(n) :=} \eqn{\delta(n)},
  where \deqn{\delta(n) =  \log\Gamma(n + 1) - n\log(n) + n - \log(2 \pi n)/2.}{%
               delta(n) =   log\Gamma(n + 1) - n*log(n) + n -  log(2*pi*n)/2.}

  Partly, pure \R transcriptions of the C code utility functions for
  \code{\link{dgamma}()}, \code{\link{dbinom}()}, \code{\link{dpois}()}, \code{\link{dt}()},
  and similar \dQuote{base} density functions by Catherine Loader.

  These \pkg{DPQ} versions typically have extra arguments with defaults
  that correspond to \R's Mathlib C code hardwired cutoffs and tolerances.
}
\usage{% >>> ../R/dgamma.R <<<
stirlerr(n, scheme = c("R3", "R4.4_0"),
         cutoffs = switch(scheme
                        , R3     = c(15, 35, 80, 500)
                        , R4.4_0 = c(5.4, 7.5, 8.5, 10.625, 12.125, 20, 26, 60, 200, 3300)
                          ),
         use.halves = missing(cutoffs),
         lgamma1p = lgamma1pC,
         order = NA,
         verbose = FALSE)

stirlerrC(n, version = c("R3", "R4..1", "R4.4_0"))

stirlerr_simpl(n, minPrec = 128L)

lgammacor(x, nalgm = 5, xbig = 2^26.5)% --> ../R/utils.R <-> ../src/lgammacor.c
}
\arguments{
  \item{x, n}{\code{\link{numeric}} (or number-alike such as "mpfr").}
  \item{verbose}{logical indicating if some information about the
    computations are to be printed.}
  \item{version}{a \code{\link{character}} string specifying the version of
    \code{stirlerrC()}.}
  \item{scheme}{a \code{\link{character}} string specifying the
    \code{cutoffs} scheme for \code{stirlerr()}.}
  \item{cutoffs}{an increasing numeric vector, required to start with
    with \code{cutoffs[1] <= 15} specifying the cutoffs to switch from 2 to 3 to ...,
    up to 10 term approximations for non-small \code{n}, where the direct
    formula loses precision.  When missing (as by default), \code{scheme}
    is used, where \code{scheme = "R3"} chooses (15, 35, 80, 500), the
    cutoffs in use in \R versions up to (and including) 4.3.z.}
  \item{use.halves}{\code{\link{logical}} indicating if the full-accuracy
    prestored values should be use when
    \eqn{2n \in \{0,1,\dots,30\}}{2n in {0,1,..,30}}, i.e.,
    \eqn{n \le 15}{n <= 15} and n is integer or integer +
    \eqn{\frac{1}{2}}{1/2}.
    Turn this off to judge the underlying approximation accuracy by
    comparison with MPFR.   However, keep the default \code{TRUE} for
    back-compatibility.}
  \item{lgamma1p}{a \code{\link{function}} used for computing
    \code{\link{lgamma}(n+1)} whenever \code{n < cutoff[1]}.  Originally
    hard wired to \code{function(n) lgamma(n + 1)} which was suboptimal for
    \eqn{0 < |n| < 1/2}.}
  \item{order}{approximation order, \code{1 <= order <= 13} or \code{NA}
    for \code{stirlerr()}.  If not \code{NA}, it specifies the number of
    terms to be used in the Stirling series which will be used for all
    \code{n}, i.e., \code{scheme}, \code{cutoffs}, \code{use.halves}, and
    \code{lgamma1p} are irrelevant.}
  \item{minPrec}{a positive integer; for \code{stirlerr_simpl} the minimal
    accuracy or precision in bits when \code{\link[Rmpfr]{mpfr}} numbers
    are used.}
  % lgammacor():
  \item{nalgm}{number of terms to use for Chebyshev polynomial approximation
    in \code{lgammacor()}.  The default, 5, is the value hard wired in \R's
    C Mathlib.}
  \item{xbig}{a large positive number; if \code{x >= xbig}, the simple
    asymptotic approximation \code{lgammacor(x) := 1/(12*x)} is used.  The
    default, \eqn{2^{26.5} = 94906265.6}, is the value hard wired in \R's C
    Mathlib.}
}
\details{
  \describe{% >> Now, most of this and more is in vignette
	    % ../vignettes/log1pmx-etc.Rnw   <<<<<<<<<
    \item{\code{stirlerr()}:}{Stirling's error, \code{stirlerr(n):=}
      \eqn{\delta(n)} has asymptotic (\eqn{n \to\infty}{n -> Inf}) expansion
      \deqn{\delta(n) = \frac 1{12 n} - \frac 1{360 n^3} + \frac 1{1260 n^5} \pm O(n^{-7}),}{%
	delta(n) = 1/(12 n) - 1/(360 n^3) + 1/(1260 n^5) +/- O(n^{-7}),}
      and this expansion is used up to remainder \eqn{O(n^{-23})} in \code{stirlerr(n)};
      different numbers of terms between different cutoffs for \eqn{n}, and
      using the direct formula for \eqn{n <= c_1}, where \eqn{c_1} is the first
      cutoff, \code{cutoff[1]}.

      Note that (new in 2024-01) \code{stirlerr(n, order = k)} will
      \emph{not} use \code{cutoffs} nor the direct formula (with
      \code{lgamma1p}), and allows \eqn{k \le 14}{k <= 14}, higher than
      with the currently maximal 10 cutoffs.
    }% stirlerr
  }% describe
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
\references{
  C. Loader (2000), see \code{\link{dbinom}}'s documentation.
  % ~/R/D/r-devel/R/src/library/stats/man/Binomial.Rd

  Our package vignette \emph{log1pmx, bd0, stirlerr - Probability Computations in R}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{dgamma}}, % ~/R/D/r-devel/R/src/library/stats/man/GammaDist.Rd
  \code{\link{dpois}}.
  High precision versions \code{\link[DPQmpfr]{stirlerrM}(n)} and % ../../DPQmpfr/man/stirlerrM.Rd
  \code{stirlerrSer(n,k)} in package \CRANpkg{DPQmpfr} (via the
  \CRANpkg{Rmpfr} and \CRANpkg{gmp} packages).
}
\examples{
n <- seq(1, 50, by=1/4)
st.n <- stirlerr(n) # now vectorized
stopifnot(identical(st.n, sapply(n, stirlerr)))
st3. <- stirlerr(n, "R3", lgamma1p = \(x) lgamma(x+1)) # previous default
st3  <- stirlerr(n, "R3", lgamma1p = lgamma1p)         # current default
## for these n, there is *NO* difference:
stopifnot(st3 == st3.)
plot(n, st.n, type = "b", log="xy", ylab = "stirlerr(n)")
st4 <- stirlerr(n, "R4.4_0", verbose = TRUE) # verbose: give info on cases
## order = k = 1:13  terms in series approx:
k <- 1:13
stirlOrd <- sapply(k, function(k) stirlerr(n, order = k))
matlines(n, stirlOrd)
matplot(n, stirlOrd - st.n, type = "b", cex=1/2, ylim = c(-1,1)/10, log = "x",
        main = "stirlerr(n, order=k) error, k = 1..13")

matplot(n, abs(stirlOrd - st.n), type = "b", cex=1/2, log = "xy",
        main = "| stirlerr(n, order=k) error |")
mtext(paste("k =", deparse(k))) ; abline(h = 2^-(53:51), lty=3, lwd=1/2)
colnames(stirlOrd) <- paste0("k=", k)

stCn <- stirlerrC(n)
all.equal(st.n, stCn, tolerance = 0)  # see 6.7447e-14
stopifnot(all.equal(st.n, stCn, tolerance = 1e-12))
stC2 <- stirlerrC(n, version = "R4..1")
stC4 <- stirlerrC(n, version = "R4.4_0")


## lgammacor(n) : only defined for n >= 10
lgcor <- lgammacor(n)
all.equal(lgcor[n >= 10], st.n[n >= 10], tolerance=0)# .. rel.diff.: 4.79e-14
stopifnot(identical(is.na(lgcor), n < 10),
          all.equal(lgcor[n >= 10],
                    st.n [n >= 10], tolerance = 1e-12))

## look at *relative* errors -- need "Rmpfr" for "truth" % Rmpfr / DPQmpfr in 'Suggests'
if(requireNamespace("Rmpfr") && requireNamespace("DPQmpfr")) {
    ## stirlerr(n) uses DPQmpfr::stirlerrM()  automagically when n is <mpfr>
    relErrV <- sfsmisc::relErrV; eaxis <- sfsmisc::eaxis
    mpfr <- Rmpfr::mpfr;     asNumeric <- Rmpfr::asNumeric
    stM <- stirlerr(mpfr(n, 512))
    relE <- asNumeric(relErrV(stM, cbind(st3, st4, stCn, stC4, lgcor, stirlOrd)))
    matplot(n, pmax(abs(relE),1e-20), type="o", cex=1/2, log="xy", ylim =c(8e-17, 0.1),
            xaxt="n", yaxt="n", main = quote(abs(relErr(stirlerr(n)))))
    eaxis(1, sub10=2); eaxis(2); abline(h = 2^-(53:51), lty=3, col=adjustcolor(1, 1/2))
    axis(1, at=15, col=NA, line=-1); abline(v=15, lty=2, col=adjustcolor(1, 1/4))
    legend("topright", legend=colnames(relE), cex = 3/4,
           col=1:6, lty=1:5, pch= c(1L:9L, 0L, letters)[seq_len(ncol(relE))])
}# Rmpfr
}
\keyword{arith}
\keyword{distribution}
