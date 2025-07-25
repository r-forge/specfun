% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/DPQ/inst/NEWS.Rd"))
\name{NEWS}
\title{DPQ News}
\encoding{UTF-8}

\section{Version 0.6-1 (2025-07-17, svn r343)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{pnchisq_ss()} gets new optional argument \code{ssr = ss(....)}
      for flexibility.

      \item New \code{gam1()} calls mpfr-ized direct code or the
      double-only C code \code{gam1d()}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{ldexp({}, 1)} or \code{pow({}, 2)} no longer segfault.
    }
  }
}

\section{Version 0.6-0 (2025-07-08, svn r342)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item added functions \code{dltgammaInc()} and \code{lgammaP11()},
      both from TOMS Algorithm 1006 by Rémy Abergel and Lionel Moisan
      (2020).

      \item C-level \code{logcf()} now uses \code{long double}, i.e., \R's
      \code{LDOUBLE}, to possibly improve for smallish \code{eps}.

      \code{logcf()} also got an optional \code{maxit} argument. For the
      pure-\R version, renamed the simplistic and slow one
      from \code{logcfR()} to \code{logcfR_vec()} and the smarter faster
      one from \code{logcfR.()} to \code{logcfR()}.  They have given
      identical results in test cases, and \code{logcfR()} should really be
      the one used practically.

      \item \code{bd0C()} has new \code{version = "R_2025_0510"} where as
      \code{bd0()} has just been updated to the latest version also in
      R-devel svn revisions 88190,1,3 and 88208.

      New \code{bd0_p1l1()}, also as \code{version} in \code{dpois_raw()};
      \code{bd0()}, and the now four \code{bd0_<meth>()} functions work
      with \code{x=Inf} and properly vectorize.

      \item The log1pmx vignette % ../vignettes/log1pmx-etc.Rnw
      now shows that for \code{dpois(x, 48)} at least, using the two parts of
      \code{ebd0()} gives clearly the most accurate values (used in R since ca. 2022).

      On the other hand, new evaluations at the end of
      \file{tests/bd0-tst.R} indicate that increasing the (only in
      \pkg{DPQ}) optional cutoff from \code{delta = 0.1} to about
      \code{delta = 0.4} gets much more accurate \code{bd0(x, np)} values, notably
      when \code{|x - np|} is smallish for relatively large \code{x} and \code{np} where
      \code{ebd0()} is clearly less accurate.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \command{rchk} \code{PROTECT} fix in \file{bd0.c}.

      \item Strictly \emph{not} comparing \code{Rboolean} with
      \code{NA_LOGICAL}.

      \item \code{log1pmx(Inf)} is correct (= \code{-Inf}) now.

      \item \code{dpois_raw(1, {})} now works (returning \code{numeric(0)}).
    }
  }
  \subsection{MISC}{
    \itemize{
      \item \file{tests/pow-tst.R} also work on "bizarre" platforms where
      \code{osVersion} is \code{NULL}.

      \item \file{tests/stirlerr-tst.R} uses many more \code{if(doExtras) ...} to run faster
      without \dQuote{extras}.
    }
  }
}

\section{Version 0.5-9 (2024-08-23, svn r333)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{dbinom_raw()} gets a new more accurate
      \code{version = "R4.4"}, using the new \code{pow1p()} function, from
      R-devel, i.e., \R 4.4.0 and newer.

      \item \code{stirlerr_simpl()} gets argument \code{version}, a string specifying
      the exact function to be used, currently with options \code{("R3",
      "lgamma1p", "MM2", "n0")}.  The first being the direct formula
      hardwired in \R up to versions 4.3.z, \code{"lgamma1p"} is to use
      \code{lgamma1p(n)} instead of \code{lgamma(n+1)}, whereas
      \code{"MM2"} is a version using \code{lgamma(n)} and slightly more
      accurate for \eqn{n \in [1, 5]}, approximately.

      \code{stirlerr()} correspondingly gets argument \code{direct.ver} to
      be passed to \code{stirlerr_simpl(., version=*)} in the
      \dQuote{direct} formula case, when \code{n < cutoffs[1]}.

      \item \code{stirlerr()} gets argument \code{order} allowing currently
      up to order 20 terms (constants \code{S0}..\code{S19}) in Stirling's
      series (and \emph{not} using \code{cutoffs} or direct formula).

      \item New \code{stirlerrC()} interfacing to C code which should
      \dQuote{parallelize} R's own C level \code{stirlerr()}.

      More documentation and tests; experiments added in
      \file{../Misc/stirlerr-trms.R}; nice explorations to \file{../tests/stirlerr-tst.R}.

      \item New functions \code{pow(x,y)} and \code{pow_di(x,y)} computing
      \eqn{x^y} as in \R's C API \file{Rmathlib}.

      \item \code{dgamma.R(x, *)} gets new optional argument \emph{dpois_r_args}.

      \item New function \code{expm1x()} providing \eqn{exp(x) - 1 - x}
      numerically accurately also for small \eqn{|x|}.

      \item New experimental \code{Ixpq()} implementing Gil et al.(2023)'s
      direct incomplete Beta function, i.e., \code{pbeta()} version.
      \bold{TODO}: make it work with \code{mpfr()} !

      \item New function \code{pntGST23_T6()} and \code{pntGST23_1()}
      implementing Gil et al.(2023)'s first asymptotic formula (Theorem 6),
      and based on Temme's Maple code for their direct computation (1).

      \code{pnt*():} CDF of the non-central t-distribution.

      \item New function \code{pntVW13()} implementing Viktor Witkovský(2013)'s
      algorithm, from his published Matlab code.

      \item Data set \code{pt_Witkovsky_Tab1} containing (the free software
      subset of) Viktor \enc{Witkovský}{Witkovsky}(2013)'s \emph{Table 1} of
      \code{pt(x, df=nu, ncp=delta)} examples; updated with octave results
      for both the original 2013 matlab code and the 2022 corrected one.

      \item \code{lgamma1p_series(x, k)}, the Taylor series approximations
      for \code{lgamma1p(x)}\eqn{ := \log \Gamma(1+x)} are now available up
      to \eqn{k = 15}, also for \code{mpfr()} numbers \code{x}.

      \item New \code{gam1d()} and \code{gamln1()} providing more TOMS 708
      auxiliary routines.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \emph{fix} format fixes for when 'long != long long'.

      \item \code{dgamma.R(x, ..)} now \emph{does} vectorize in \code{x},
      as documented.

      \item minor adjustment of \code{pntR()}'s default \code{use.pnorm}.

      \item \file{src/bd0.c} needed one more \code{PROTECT()}.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item new \file{tests/dt-ex.R} exploring \code{dt()} accuracy issues,
      mostly historical, but newly for \code{df=}\eqn{\nu \ll 1}{nu << 1}.
    }
  }
}

\section{Version 0.5-8 (2023-11-30, svn r278)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{ebd0()} and C level \code{ebd0C()} now work with \emph{long} vectors and
      hence return \code{data.frame(yh = *, yl = *)} instead of 2-row matrix.

      \item format fixes: \code{"\%lld"} and \code{(long long)} notably for
      \code{R_xlen_t}.
    }
  }
}

\section{Version 0.5-7 (2023-10-20, svn r275)}{
  \subsection{Misc}{
    \itemize{
      \item fix "lost braces"
      \item replace UTF-8 greek glyps by plotmath
    }
  }
}

\section{Version 0.5-6 (2023-09-06, svn r272)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{phyperR()} gets \code{lower.tail=TRUE, log.p=FALSE} arguments.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Correctly vectorize \code{phyperR()}.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item Documentation (and internal comment) updates for
      \code{qbinomR()} and \code{qnbinomR()}, notably fixing "empty"
      argument doc.

      \item Deprecated \code{phyperBinMolenaar()}; it has been
      synonym to \code{phyperBinMolenaar.1()}.

      \item Fix \LaTeX warnings/errors by removing amazingly many
      \verb{\cr} in \file{man/*.Rd} files.
    }
  }
}

\section{Version 0.5-5 (2023-06-23, svn r264)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{qbinomR(p, *)}, \code{qnbinomR(p, *)}, and \code{qpoisR(p, *)}:
      Fix \code{R_Q_P01_boundaries(p, ..)} for trivial boundaries, e.g.,
      \eqn{p \in \{0,1\}} for \code{log.p=FALSE}.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item A new \command{flang} Fortran compiler does not know outdated
      \code{derfc()}, but needs current standard \code{erfc()}.
    }
  }
}

\section{Version 0.5-4 (2023-04-12, svn r261)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \R function \code{qnormCappr()}.
      \item New \R function \code{rexpm1()} from TOMS 708.
      \item Provide also \code{qnormR(*, version = "1.0.x")} and \code{"1.0_noN"},
      for history's sake.
      \item New \R function \code{rlog1()} from TOMS 708.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{qnormUappr()} & \code{qnormUappr6()} were \emph{wrongly}
      negative for lp in \eqn{[-.693147, 0)} (corresponding to \eqn{p > 1/2}).
    }
  }
  \subsection{Misc}{
    \itemize{
      \item Tweaks to the \file{qnorm-asymp} vignette, related to the JSS submission.

      \item \file{tests/pqnorm_extreme.R} plots more about the
      \code{pnormAsymp()} relative errors.

      \item Apple clang 14.0.3 needs looser tolerance for \code{logcf{R}()}
      and \code{dbinom_raw()} comparisons.

      \item Minimally document the \emph{internal} \code{.D_*()} and
      \code{.DT_*()} utilities, as we provide them to experts and e.g.,
      package \CRANpkg{DPQmpfr}.
    }
  }
}

\section{Version 0.5-3 (2022-12-01, svn r240)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New \R functions \code{qntR()} and \code{qtU()},
      \code{Vectorize()}d versions of \code{qntR1()} and \code{qtU1()},
      respectively.

      \item New \code{qtR() <- Vectorize(qntR1, *)}; pure \R
      implementation of \R's Mathlib C-level \code{qt()}, but additionally
      allowing tweaks, used for fixing \R's bug \PR{18630}.
      Added (optional) log-scale Newton steps via \code{logNewton = log.p},
      needing more tests.

      \item New \code{qtNappr()} -- from the remark in R's \file{nmath/qt.c}
      about very large \code{df} approximation, now with all 4 large-df terms
      from Abramowitz & Stegun's (26.7.5).

      \item New \code{gammaVer()} to exemplify \R's (partly historical) % ../man/gammaVer.Rd
      versions of \code{gamma()} implementations.

      \item New \code{qnormUappr6()}, % ../man/qnormUappr.Rd
      providing the \sQuote{6 coefficients}-approximation of A. & S. to \code{qnorm()}.

      \item \code{qnormR()} gets new \code{version = "2022-08-04"} which
      uses MM's newly derived asymptotic approximations to \code{qnorm()}.

      \item New \code{qnormAsymp()} providing new asymptotic approximations
      to \code{qnorm()}.

      \item Added vignette \dQuote{Asymptotic Tail Formulas for Gaussian Quantiles}.

      \item New \code{chebychevPoly()} (and auxiliaries) for evaluation
      of Chebyshev polynomials.
    }
  }
}

\section{Version 0.5-2 (2022-06-08, svn r212)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{qbetaAppr()} (and the \code{qbetaAppr.{1,3}()}
      auxiliaries) now get a \code{log.p} argument.

      \item all(?) \code{qbeta.Appr*()} functions now accept a
      \code{lower.tail} argument (with default \code{TRUE}). % __FIXME__ Needs TESTS!
      %> ../tests/qbeta-dist.R -- C-s for 'FIXME: More qbetaAppr*' (and
      %%%lower.tail ..)

      \item New \R function \code{dpsifn()}, as interface to \R's C API
      \code{dpsifn()}, the workhorse of all \R's \code{psigamma()}
      functions.

      \item \code{dpois_raw(x, lambda, *)} gets new cutoff \code{small.x__lambda}
      to \emph{not} use \code{bd0()} or \code{ebd0()} when \code{x} is much
      smaller than \code{lambda}.

      \item Simple functions, mostly for didactical and comparison purposes,
      \code{dpois_simpl()}, \code{dpois_simpl0()} and
      \code{stirlerr_simpl()}, all of which use the \dQuote{simple} direct
      formula which leads to numerical cancellation typically.

      \item New \R function \code{bpser()} interfacing to the C function of
      the same in \R's Mathlib \code{pbeta()} or \code{bratio()} in file
      \file{nmath/toms708.c}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Updated URL of Abramowitz & Stegun.
      \item Fixed \code{\\} escape in one \file{man/*.Rd}
    }
  }
  \subsection{Misc}{
    \itemize{
      \item Our \file{tests/*.R} no longer rely on the \CRANpkg{Matrix}
      \file{test-tools.R} collection but rather use our own. % ./test-tools.R

      \item In \file{tests/chisq-nonc-ex.R} (and \file{DESCRIPTION}), we've
      replaced \pkg{akima} by \CRANpkg{interp} as the latter has a FOSS licence.
    }
  }
}

\section{Version 0.5-1 (2021-12-10, svn r197)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{bd0(x, np, delta, ..)} now uses
      \code{if (|x-np| <= delta * (x+np))}  (less-equal instead of strictly
      less), and hence setting \code{delta = 0} is now allowed, using the
      series expansion only for \code{x == np}, useful e.g., for the case
      of highly accurate mpfr-numbers.

      Similarly, in \code{log1pmx(x, .., eps2, ..)},  \code{eps2 = 0} is now allowed.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item consistency between \code{ebd0()} and \code{ebd0C()}:
      multiplication \code{ * e } coming last ("bug" did not show ..).

      \item \code{okLongDouble()} should no longer fail on M1 mac.

      \item \code{ebd0C()} now checking \eqn{|yl| < 5.5} before returning; was
      horrendously wrong on Windows compiled with \command{-mnative} as long
      as it included (the default) \command{-mfma} (FMA:= Fused-Multiply-Add).

      \item fix \command{rchk} issue: "need" \code{PROTECT(.)} while calling allocating \code{lgamma1p()}.
    }
  }
}

\section{Version 0.5-0 (2021-09-10, svn r183)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item new \code{ebd0C()} interface to C version; helps to fix
      lapsus in pure \R version \code{ebd0()}:

      \code{ebd0(x, M)} for large \code{M} now checks for overflow in
      \code{M/x}, and notably the case \code{f * 2^-(e+10) =: fg == Inf}.

      \item In \code{bd0()} work around underflow of \code{(x-np)/(x+np)}.

      \item new \R functions \code{frexp()} and \code{ldexp()} for getting
      and setting base-2 representations of numbers,
      and new \R function \code{modf()} to split number into integer and fractional part,
      all interfacing the C99 (math lib) standard functions of the same name.

      \item \code{logcfR.(x)} has been \dQuote{vectorized} in \code{x} even
      though it's iterative with different number of iterations for each
      \code{x[i]}, and is hence considerably faster when \code{x} is an
      \code{"mpfr"} vector.

      TODO: Consider renaming the two logcfR versions; at least neither is deprecated!
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item fixed embarrassing typo (two \code{i}'s) bug in \code{logcfR()}

      \item \code{stirlerr(n)} now also works when \code{is.integer(n)} and
      \code{n} is large enough for \code{n * n} to overflow (to \code{NA}).

      \item fix typo/thinko in \code{dnbinomR()}

      \item ditto in \code{logcfR(*, trace=TRUE)} iteration report.

      \item \code{dpois_raw(x, *)} now works up to maximal \code{x}, now
      preventing overflow in previous \code{2*pi*x} computation.

      \item \code{dnbinomR()} and \code{dnbinom.mu()} fix for \code{x >
	size} when \code{log=TRUE}, notably for \code{x >> size}.

    }
  }
}

\section{Version 0.4-4 (2021-05-22, svn r175)}{% CRAN publication 05-27
  \subsection{NEW FEATURES}{
    \itemize{
      \item new \code{bd0_*()} versions of \code{bd0()}, based mainly on \code{log1pmx()}.
    }
  }
  \subsection{REFACTORING}{
    \itemize{
      \item To be more modular, our \file{test-tools.R} no longer
      \code{source()} those of \pkg{Matrix}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \file{src/bd0.c} gave warning with some compiler settings with
      \code{-Wself-assign}.% but our gcc 10.2.1 (20201125) does not work with -Wself-assign

      \item \file{test-tools.R} \code{readRDS_()} thinko fixed.
    }
  }
}

\section{Version 0.4-3 (2021-05-05, svn r166)}{% = CRAN publication
  \subsection{NEW FEATURES}{
    \itemize{
      \item Provide \R functions
      \code{log1mexpC()}, \code{log1pexpC()}, \code{log1pmxC()}, and \code{lgamma1pC()},
      all interfacing to \R's C API (\file{Rmath.h}), aka \sQuote{Rmathlib}.

      \item New \code{pnormAsymp()} for asymptotic (typically upper tail,
      log scale) approximation of \code{pnorm()}.  With \CRANpkg{Rmpfr}, we
      can see how accurate these are in detail. % currently on  ../man/pnormAsymp.Rd

      \item New \code{dnbinomR()} finding better code for \R itself,
      \code{dnbinom.mu()} and \code{dbinom_raw()}; also new utility
      functions \code{bd0()} and \code{stirlerr()}, both vectorized,
      and also working with "mpfr"-numbers, such that \code{dpois_raw()}
      now does so, too.

      Additionally, an \dQuote{extended} version of \code{bd0()} called
      \code{ebd0()} in pure \R, where a C version was proposed by Morten
      Welinder in R's bugzilla, \PR{15628}. % However, ebd0() is *wrong* too

      Experimentally, also provide \code{p1l1()} and its Taylor series
      approximations \code{p1l1ser()} which could be employed for
      \code{bd0()} / \code{ebd0()} instead of the current algorithms.

      \item Several C level utilities to be .Call()ed from \R, from \R's
      mathlib, \code{ldexp()} and \code{frexp()} even from C math.

      \item Pure R implementations (somewhat experimental) of corresponding
      \R Mathlib C functions: \code{qbinomR()}, \code{qnbinomR()}, and \code{qpoisR()}
      each with several tuning parameters for the underlying algorithm,
      notably the root-finding parts.

      \item newly, \code{logcf()} now based on C code, perfectly
      vectorizes; the pure \R version, now called \code{logcfR(x, *)}
      currently still runs the iterations simultaneously for all `x[i]` and
      hence convergence and rescaling happen by \dQuote{group decision},
      efficiently but undesirable for strict comparisons.

      \code{logcfR(x, *)} and \code{log1pmx(x, *)} now also work for
      "mpfr"-numbers \code{x}, and \code{log1pmx()} gets optional
      \code{eps2 = 1e-2} and \code{minL1 = -0.791} arguments, the defaults
      of which \emph{may change}, as I think at least the latter to not be
      perfect.

      \item Now \code{lb_chiAsymp(nu, order)} works up to order 8.

      \item Provide the first parts of a new vignette
      (\file{../vignettes/log1pmx-etc.Rnw}) on \code{log1pmx()},
      \code{bd0()}, and \code{stirlerr()}, which should contain part of
      Loader(2000)'s report and new findings of improved \code{bd0()} and
      \code{stirlerr()} computations.
    }
  }
}

\section{Version 0.4-2 (2020-11-07, svn r151)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New functions \code{pnormL*()} and \code{pnormU*()} for
      (mathematically proven) lower and upper bound to \code{pnorm()},
      notably also for investigation with \code{log.p=TRUE}.

      \item \code{qnormR()}, implementing current \R's \code{qnorm()} in
      pure R, notably with \code{trace} and \code{version} options.
    }
  }
}

\section{Version 0.4-1 (2020-06-17, svn r145)}{
  \subsection{TESTING}{
    \itemize{
      \item Reorganize tests; notably to become less platform dependent.
    }
  }
}

\section{Version 0.4-0 (2020-06-15, svn r137)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Many new \code{phyper*()} functions and helpers for them, such
      as Bernoulli numbers \code{Bern()} and asymptotic
      \code{lgammaAsymp()}.

      \item Notably \code{phyperR2()} which is a pure \R version of \R's
      own (C code based) \code{phyper()}.
    }
  }
}

\section{Version 0.3-5 (2019-10-18, svn r131)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{pnbetaAS310()} gained a \code{useAS226} option.

      \item New \code{okLongDouble()} function, notably for detecting that
      with a \command{valgrind}ed version of \R, \code{long double} C
      arithmetic seems to silently switch to (simple) \code{double} precision.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item long double printing from C now should happen correctly on all
      platforms including 32- and 64-bit Windows.
    }
  }
}

\section{Version 0.3-4 (2019-10-16, svn r125)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item more efficient (internally vectorized) \code{dntJKBf()}.
      Consequently, \code{dntJKBf1()} is deprecated now.

      \item \code{pntR()} (and \code{pntR1()}) get new option
      \code{use.pnorm} (the default of which had been hard coded previously).
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item fix thinko in \code{any_mpfr()} and \code{all_mpfr()}.

      \item \code{pnchisqRC()}'s C code gave severe valgrind warnings;
      fixed printing of \code{long double} etc;
      also added special MinGW deal in Windows.

      \item \code{ppoisD()} behaves differently in a valgrinded version;
      for now, reproducible only when using valgrind on non-valgrinded
      installed package.
    }
  }
}


\section{Version 0.3-3 (2019-09-24, svn r110)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Renamed (and generalized / merged) many functions to have less "." in names.

      \item New \code{pnbetaAS310()} function interfacing to my
      \emph{corrected} C version of 'ASA 310' (2007).

      \item New \code{algdiv()} function interfacing to the 'TOMS 708' version of
      our \code{logQab()}.

      \item New \code{pl2curves()} which generalizes (somewhat) previous
      function \code{p.dnchiB()}.

      \item Made \code{newton()} more flexible with new \code{xMin} and
      \code{xMax} arguments which are notably useful for \code{q*()}
      (quantile function) computations.

      Correspondingly replaced previous \code{qchisq2()} and
      \code{qchisqA()} with new \code{qchisqN()} (\sQuote{N}ewton).

      \item new \code{pnchisqRC()}, a version of R's C level non-central
      chi-squared, with additional options.

      \item \dQuote{new} \code{logspace.add()} and \code{logspace.sub()}.

      \item \dQuote{new} \code{pnchisqT93()} (plus its two auxiliaries),
      implementing Temme(1993)'s approximations.

      \item \dQuote{new} \code{pnchisqBolKuz()} and \code{qnchisqBolKuz()}
      implementing Bol'shev and Kuznetzov (1963)'s approximations.

      \item \dQuote{new} \code{pnchi1sq()} and \code{pnchi2sq()} with
      \dQuote{exact} formulas for the special cases \code{df=1} and
      \code{df=3}.

      \item simplified formula in \code{dtWV()}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{qnchisqPearson(pp, df=DF, ncp=100)} and
      \code{qnchisqSankaran_d(*)} no longer return \code{NaN} for very
      large \code{DF = 1e200}.

      \item \code{pnchisq()} now also has default \code{verbose = 0} as all
      other such functions with \code{verbose} (or \code{trace[.lev]}
      optional argument.
    }
  }
}

\section{Version 0.3-0 [2018-08-28]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Move many of my up to 15 years old DPQ computation utilities
      into a package, to become public, \dQuote{bloggable}, etc.
    }
  }
}
