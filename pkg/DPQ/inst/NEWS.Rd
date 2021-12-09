% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/DPQ/inst/NEWS.Rd"))
\name{NEWS}
\title{DPQ News}
\encoding{UTF-8}

\section{Version 0.5-1 (2021-12-09, svn r196)}{
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
